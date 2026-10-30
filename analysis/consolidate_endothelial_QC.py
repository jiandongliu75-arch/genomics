"""Consolidate strict and explicitly separate nuclear-supported EC candidates."""
import h5py
import anndata as ad
from anndata.io import read_elem, write_elem
from endothelial_common import *
from reassess_endothelial_RNA import collapse_read

def refresh_gene_names():
    """Exclude technical gene aliases without changing a graph, family or count."""
    folder=OUT/'snRNA'
    rows=pd.read_csv(folder/'EC_cluster_summary.tsv',sep='\t')
    markers=pd.read_csv(folder/'EC_cluster_markers.tsv',sep='\t')
    technical=markers.gene.str.lower().str.startswith('mt-') | markers.gene.str.match(r'^Rp[sl]\d')
    markers=markers[~technical].copy()
    for i,r in rows.iterrows():
        d=markers[markers.cluster_id.eq(r.cluster_id)]
        both=d[d.positive_libraries.eq(2)].gene.tolist()
        top=both if both else d.gene.tolist()[:3]
        rows.loc[i,'gene_name']='/'.join(top[:2])+' ECs' if top else 'Marker-poor ECs'
        rows.loc[i,'top_markers']=';'.join(top[:12]);rows.loc[i,'two_library_markers']=len(both)
        keys=['QC_adjusted_cluster_best_jaccard','markers_withheld_cluster_best_jaccard','seed_7_best_jaccard','seed_123_best_jaccard']
        ok=r.nuclei>=50 and r.libraries_ge10==2 and r.largest_library_fraction<=.9 and len(both)>=3 and r[keys].min()>=.5
        rows.loc[i,'grouping_support']='Supported RNA grouping' if ok else 'Candidate / sensitivity-dependent grouping'
    names=rows.set_index('cluster_id').gene_name;support=rows.set_index('cluster_id').grouping_support
    save(rows,folder/'EC_cluster_summary.tsv');save(markers,folder/'EC_cluster_markers.tsv')
    ann=pd.read_csv(folder/'EC_annotations.tsv.gz',sep='\t',index_col=0)
    ann['EC_gene_name']=ann.EC_cluster_id.map(names);ann['EC_grouping_support']=ann.EC_cluster_id.map(support)
    ann.to_csv(folder/'EC_annotations.tsv.gz',sep='\t',na_rep='NA')
    for stage in STAGES:
        a=ann[ann.age.eq(stage)]
        a.to_csv(folder/stage/'EC_annotations.tsv.gz',sep='\t',na_rep='NA')
        for file in ['EC_corrected_and_raw_counts.h5ad','EC_embeddings_metadata_only.h5ad']:
            with h5py.File(folder/stage/file,'r+') as h:
                old=read_elem(h['obs']);np.testing.assert_array_equal(old.index,a.index)
                old['EC_gene_name']=a.EC_gene_name;old['EC_grouping_support']=a.EC_grouping_support
                del h['obs'];write_elem(h,'obs',old)
    p=pd.read_csv(folder/'EC_pseudobulk_metadata.tsv',sep='\t',index_col=0)
    p['EC_gene_name']=p.EC_cluster_id.map(names);p.to_csv(folder/'EC_pseudobulk_metadata.tsv',sep='\t')

def threshold_sensitivity(obs,x,genes,qc):
    ev=lineage_evidence(x,genes,'snRNA');ev.index=obs.index;rows=[]
    coherent=qc.coherent_mural_group_review.to_numpy()
    for age in STAGES:
        m=obs.age.eq(age).to_numpy()
        for min_ec in [1,2,3]:
            for foreign in [.01,.02,.03]:
                ec=(ev.Endothelial_genes>=min_ec)&(ev.Endothelial_fraction>=.001)
                mixed=np.zeros(len(obs),bool)
                for c in LINEAGES:
                    if c!='Endothelial':mixed|=(ev[c+'_genes']>=3)&(ev[c+'_fraction']>=foreign)&(ev[c+'_molecules']>=10)
                rows.append({'stage':age,'minimum_EC_markers':min_ec,'foreign_fraction_threshold':foreign,
                    'clean_before_mural_group_review':int(np.sum(m&ec&~mixed)),
                    'primary_clean_nuclei':int(np.sum(m&ec&~mixed&~coherent)),
                    'weak_EC_review':int(np.sum(m&~ec)),'mixed_program_review':int(np.sum(m&ec&mixed)),
                    'coherent_mural_group_review':int(np.sum(m&ec&~mixed&coherent)),
                    'mural_QC_groups_held_fixed':True,'primary_rule':min_ec==2 and foreign==.02})
    return pd.DataFrame(rows)

def main():
    refresh_gene_names()
    obs,var,genes,x,raw=collapse_read()
    qc=pd.read_csv(OUT/'snRNA/all_endothelial_QC_decisions.tsv.gz',sep='\t',index_col=0)
    ann=pd.read_csv(OUT/'snRNA/EC_annotations.tsv.gz',sep='\t',index_col=0)
    sensitivity=threshold_sensitivity(obs,x,genes,qc);marker_rows=[]
    # Explicit sensitivity around the conservative residual-program screen.
    ev=lineage_evidence(x,genes,'snRNA');ev.index=obs.index
    for age in STAGES:
        m=obs.age.eq(age).to_numpy()
        ix=np.flatnonzero(m & ~qc.primary_EC_eligible.to_numpy())
        review=ad.AnnData(x[ix],obs=qc.iloc[ix].copy(),var=var.copy());review.layers['counts_raw']=raw[ix]
        review.uns['purpose']='QC review nuclei retained for audit, excluded from primary endothelial subclustering. Mixed RNA programs are not proven doublets.'
        review.write_h5ad(OUT/'snRNA'/age/'EC_review_corrected_and_raw_counts.h5ad',compression='gzip')
    for cid,group in ann.groupby('EC_cluster_id',observed=True):
        for lib,sub in group.groupby('sample',observed=True):
            ii=obs.index.get_indexer(sub.index)
            for layer,counts in [('DecontX_corrected',x),('raw',raw)]:
                z=lognorm(counts[ii]);f=np.asarray((counts[ii]>=1).mean(axis=0)).ravel();mean=np.asarray(z.mean(axis=0)).ravel()
                for fam,panel in PANELS.items():
                    ix=np.flatnonzero(np.isin(genes,panel));detected=genes[ix[(f[ix]>=.15)&(mean[ix]>=.1)]]
                    marker_rows.append({'stage':str(sub.age.iloc[0]),'cluster_id':cid,'sample':lib,'layer':layer,'nuclei':len(ii),
                        'family':fam,'detected_markers':';'.join(detected),'detected_anchors':';'.join(sorted(set(detected)&set(ANCHORS[fam]))),
                        'detected_marker_count':len(detected),'detected_anchor_count':len(set(detected)&set(ANCHORS[fam]))})
    save(sensitivity,OUT/'snRNA/QC_threshold_sensitivity.tsv')
    save(pd.DataFrame(marker_rows),OUT/'snRNA/family_markers_by_library.tsv')
    save(ann.groupby(['age','sample','EC_family'],observed=True).size().reset_index(name='nuclei'),OUT/'snRNA/EC_family_counts_by_library.tsv')
    combined=[];summaries=[];cohort_rows=[]
    for stage in STAGES:
        folder=OUT/'Xenium'/stage
        a=pd.read_csv(folder/'EC_consensus_annotations.tsv.gz',sep='\t',index_col=0,float_precision='round_trip')
        nuclei=pd.read_parquet(BASE/'results/xenium_boundary_inputs'/stage/'nucleus_boundaries.parquet',columns=['cell_id','label_id']).drop_duplicates()
        components=nuclei.groupby('cell_id',observed=True).size().reindex(a.index).fillna(0).astype(int)
        geometry=pd.read_csv(ROOT/'results/xenium_nuclear_anchored_segmentation'/stage/'segmentation_profile_audit.tsv.gz',sep='\t',index_col=0).reindex(a.index)
        a['deposited_projected_nuclear_components']=components
        a['multiple_projected_nuclear_components_review']=components.gt(1)
        a['segmentation_seed_type']=geometry.seed_type
        a['segmented_area_um2']=geometry.segmented_area_um2
        a['native_to_refined_area_ratio']=geometry.segmented_area_um2/a.cell_area.replace(0,np.nan)
        a['native_nuclear_area_missing']=a.nucleus_area.isna()
        merged=a.multiple_projected_nuclear_components_review
        a.loc[merged,'confidence_tier']='Multiple projected nuclear components; merge review'
        a.loc[merged,'fine_consensus_family']=UNRESOLVED
        a.loc[merged,['high_confidence_subtype','assigned_including_moderate','EC_broad_clean']]=False
        a['new_broad_EC_candidate']=~a.native_EC_broad & a.EC_broad_clean
        frames={v:pd.read_csv(folder/f'EC_{v}_assignments.tsv.gz',sep='\t',index_col=0) for v in
                ['native','nuclear_anchored','high_confidence','depth_matched','nuclear_only','nuclear_calibrated','nuclear_depth_matched','cytoplasmic_only']}
        nuc,nucal=frames['nuclear_only'],frames['nuclear_calibrated']
        n_ok=nuc.fine_supported & nucal.fine_supported & nuc.mapped_family.eq(nucal.mapped_family)
        # A supported incompatible whole-cell/core/native subtype blocks rescue.
        conflict=np.zeros(len(a),bool)
        for v in ['native','nuclear_anchored','high_confidence']:
            conflict|=(frames[v].fine_supported & frames[v].mapped_family.ne(nucal.mapped_family)).to_numpy()
        n_ok &= ~conflict & ~a.geometry_unresolved_review & ~merged
        a['nuclear_supported_family']=np.where(n_ok,nucal.mapped_family,UNRESOLVED)
        a['nuclear_supported_subtype']=n_ok
        a['nuclear_only_subtype_candidate']=n_ok & ~a.assigned_including_moderate
        a['nuclear_rescue_with_mixed_cell_RNA']=a.nuclear_only_subtype_candidate & frames['nuclear_anchored'].mixed_lineage_review
        a['nuclear_broad_EC_clean']=nucal.broad_EC_competition_pass & ~nucal.mixed_lineage_review & ~a.geometry_unresolved_review & ~merged
        a['stage']=stage
        # Refit the independent spatial graph after geometry QC, on the same
        # RNA-information-eligible cohort in all compared ownership views.
        from reassess_endothelial_Xenium import cluster_views
        counts=ad.read_h5ad(folder/'EC_candidate_count_views.h5ad')
        values={'nuclear_anchored':counts.X.toarray().astype(np.float32)}
        for v in ['high_confidence','nuclear_only','nuclear_depth_matched']:
            values[v]=counts.layers[v].toarray().astype(np.float32)
        common=~merged & ~a.geometry_unresolved_review & ~frames['nuclear_calibrated'].mixed_lineage_review
        for v in ['nuclear_anchored','high_confidence','nuclear_only','nuclear_depth_matched']:
            common &= frames[v].information_QC
        clusters=cluster_views(stage,a.index,values,frames,counts.var_names.to_numpy(dtype=str),folder,extra_eligible=common)
        for c in ['Xenium_cluster_id','Xenium_gene_name','Xenium_UMAP1','Xenium_UMAP2']:
            a[c]=clusters[c].reindex(a.index)
        del counts,values
        a.to_csv(folder/'EC_consensus_annotations.tsv.gz',sep='\t',na_rep='NA')
        # Refresh obs only. Keep every count view bit-for-bit intact.
        with h5py.File(folder/'EC_candidate_count_views.h5ad','r+') as h:
            old=read_elem(h['obs']);np.testing.assert_array_equal(old.index,a.index)
            del h['obs'];write_elem(h,'obs',a)
        combined.append(a)
        for cohort,m in [('all_EC_candidates',np.ones(len(a),bool)),('native_broad_EC',a.native_EC_broad.to_numpy())]:
            d=a.loc[m]
            summaries.append({'stage':stage,'cohort':cohort,'profiles':len(d),'new_cell_broad_clean':d.EC_broad_clean.sum(),
                'nuclear_broad_clean':d.nuclear_broad_EC_clean.sum(),'high_confidence_A_B':d.high_confidence_subtype.sum(),
                'assigned_A_B_C':d.assigned_including_moderate.sum(),'nuclear_supported_subtypes_all':d.nuclear_supported_subtype.sum(),
                'additional_nuclear_only_candidates':d.nuclear_only_subtype_candidate.sum(),
                'additional_nuclear_candidates_with_mixed_cell_RNA':d.nuclear_rescue_with_mixed_cell_RNA.sum(),
                'multiple_projected_nuclear_components_review':d.multiple_projected_nuclear_components_review.sum(),
                'no_projected_nuclear_component_record':d.deposited_projected_nuclear_components.eq(0).sum(),
                'strong_CM_program_native':frames['native'].loc[m,'strong_Cardiomyocyte_program'].sum(),
                'strong_CM_program_nuclear':nucal.loc[m,'strong_Cardiomyocyte_program'].sum(),
                'mixed_new_cell_review':frames['nuclear_anchored'].loc[m,'mixed_lineage_review'].sum()})
            for v,f in frames.items():
                # Native-broad EC is a frozen selection, so all views use identical IDs.
                e=f.loc[m]
                cohort_rows.append({'stage':stage,'cohort':cohort,'variant':v,'profiles':len(e),
                    'mean_CM_fraction':e.Cardiomyocyte_fraction.mean(),'median_CM_fraction':e.Cardiomyocyte_fraction.median(),
                    'mean_EC_fraction':e.Endothelial_fraction.mean(),'strong_CM_fraction':e.strong_Cardiomyocyte_program.mean(),
                    'median_molecules':e.total_gene_molecules.median(),'broad_EC_competition_pass':e.broad_EC_competition_pass.sum(),
                    'fine_supported':e.fine_supported.sum()})
        done=json.loads((folder/'complete.json').read_text())
        done['nuclear_supported_subtypes']=int(n_ok.sum());done['additional_nuclear_only_subtype_candidates']=int(a.nuclear_only_subtype_candidate.sum())
        done['high_confidence_subtypes']=int(a.high_confidence_subtype.sum());done['assigned_including_moderate']=int(a.assigned_including_moderate.sum())
        done['broad_clean_after_segmentation']=int(a.EC_broad_clean.sum());done['new_broad_candidates']=int(a.new_broad_EC_candidate.sum())
        done['multiple_projected_nuclear_components_review']=int(merged.sum())
        done['de_novo_final_graph_cohort']='Whole-cell broad EC competition pass, no strong whole-cell foreign program, no unresolved/multiple-component geometry, and information-eligible in cell/core/nuclear/nuclear-depth-control views. Final graphs refitted after these QC exclusions.'
        done['multi_nuclear_review_limit']='Multiple deposited projected nuclear components can indicate merged profiles or genuine multinucleation. Such parents are retained for review and excluded from strict fine release; no unsupported claim of successful cell splitting.'
        done['broad_assay_calibration']='Native whole-section profiles calibrate broad reference competition; EC candidate cohort calibrates fine subtypes. Nuclear sensitivity separately calibrates all-section nuclear broad features and informative candidate nuclear fine features.'
        done['nuclear_candidate_limit']='Extra nuclear-only candidates are separately reported; mixed cell RNA does not become a fully clean whole-cell assignment. Both nuclear scalings must support the same family and no supported cell/core/native family may conflict.'
        write_json(folder/'complete.json',done)
    pd.concat(combined).to_csv(OUT/'Xenium/all_stage_EC_consensus.tsv.gz',sep='\t',na_rep='NA')
    save(pd.DataFrame(summaries),OUT/'Xenium/EC_QC_summary.tsv')
    save(pd.DataFrame(cohort_rows),OUT/'Xenium/EC_RNA_ownership_depth_comparison.tsv')
    write_json(OUT/'methods.json',{'snRNA_source':str(RNA_INPUT),'snRNA_existing_filter':'scDblFinder singlets; >=1000 corrected UMIs;750-7500 genes;<=5% mitochondrial;DecontX contamination<=50%;eight converged per-library fits',
        'snRNA_EC_screen':'Whole-genome denominator: >=2 EC diagnostic genes with corrected count>=1 and EC fraction>=0.1%; review competing program if>=3 genes,>=10 molecules,>=2% full-transcriptome UMIs. All reviewed nuclei/counts retained.',
        'snRNA_two_pass_mural_review':'Provisional within-age RNA clusters with >=3 mural markers detected in>=15% of nuclei in both raw/corrected counts and both libraries (>=10 nuclei/library) are retained as coherent mixed-group review; final primary graphs/refits exclude these groups. This is not proof every nucleus is a doublet.',
        'spatial_EC_screen':'403-target denominator: >=3 EC genes and >=3% EC molecules; >=2 EC genes for nuclear views. Foreign-program review:>=3 genes,>=10 molecules,>=10% target molecules.',
        'thresholds_not_comparable_across_assays':True,'RNA_new_clustering':'Separate ages;within-age two-library Harmony;2000HVG minus mt/ribo/sex and foreign-lineage diagnostics;30PC,20NN,Leiden.45,seed42; alternative seeds,resolutions,raw counts,QC covariates,explicit annotation markers withheld,ambient<=30% reclustering.',
        'RNA_family_assignment':'Post-hoc curated and paper panels;>=3 expressed markers,>=2 anchors (capillary>=1);score>0.05,top-v-other margin>=.10;lymphatic requires Prox1 and Flt4. Unqualified top competing panels can block a label.',
        'RNA_gene_names':'Top two within-stage enriched genes, preferentially consistent direction in both libraries; aliases describe enrichment, not mandatory per-cell coexpression or verified function.',
        'Xenium_segmentation':'Reuse completed image-only whole-section DAPI-derived nuclear anchors + acquired ATP1A1/CD45/E-Cadherin boundary channel at0.425um. Native-profile territories retained,2um mask expansion,0.85um uncertain boundary band excluded in core view. Source nuclear molecule ownership conserved. No new nuclear detection or cell-count/merge-split validation.',
        'Xenium_scope':'All native EC profiles plus original whole-reference EC proposals among uncertain profiles; other cardiac profiles compete as RNA reference alternatives. New broad candidates remain labeled candidates; parent atlas not overwritten.',
        'Xenium_background':'QV>=20 gene molecules; unassigned RNA not forced into primary cells; 403 genes exclude negative-control features. No claim of global ambient removal; compare native/new/core/nuclear/cytoplasm and depth-matched controls.',
        'Xenium_QC_tiers':'A/B: new whole-cell and both nuclear scalings same supported subtype;A also native. C: native+new+boundary-core same supported subtype. Mixed whole/nuclear programs and supported subtype conflicts excluded from A/B/C. Extra nuclear-only candidates separately named and never promoted to clean whole-cell evidence.',
        'Xenium_nuclear_geometry_QC':'Count deposited (cell_id,label_id) nuclear components. Parents with>1 projected component are retained as possible merge/multinucleation review and excluded from A/B/C and nuclear-only candidate release. Missing nuclear area/outline alone never deletes a profile. Refined/native area changes and fallback seeds are exported as diagnostics, not automatically called errors.',
        'confidence_limit':'Operational robustness criteria, not calibrated probabilities or biological/anatomical accuracy; count views share molecules/references.',
        'biological_limits':['One Xenium section per age; two RNA libraries per age are not assumed independent animals.','No absence claim from zero fine calls.','CM-associated endothelial RNA can reflect genuine biology, ambient RNA, unresolved merges or ownership error; present evidence does not prove its origin.','Nuclear RNA reduces exposure to boundary-associated mixing but may contain ambient/cross-projection errors and omit real cytoplasmic state genes.','No coordinate smoothing, subtype quotas or age regression; no new causal signaling or lineage-conversion inference.'],
        'panels':PANELS,'anchors':ANCHORS,'lineage_diagnostics':LINEAGES,'segmentation_reference':'../xenium_nuclear_anchored_segmentation/segmentation_design.json'})
    print(pd.DataFrame(summaries).to_string(index=False),flush=True)

if __name__=='__main__':main()
