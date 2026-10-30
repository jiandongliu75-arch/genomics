"""Data-level validation of endothelial QC, retained counts and released labels."""
import h5py
import anndata as ad
from anndata.io import read_elem, sparse_dataset
from PIL import Image
import subprocess
import re
from endothelial_common import *
from reassess_endothelial_RNA import collapse_read

def exact_sparse(a,b):
    assert a.shape==b.shape
    diff=(a-b).tocsr();diff.eliminate_zeros()
    assert diff.nnz==0, (diff.nnz,float(np.max(np.abs(diff.data))))

def compact_h5ad(path):
    """Remove HDF5 free-space left by obs-only updates; preserve stored arrays."""
    before=path.stat().st_size
    temp=path.with_name(path.name+'.repack')
    try:
        with h5py.File(path,'r') as source,h5py.File(temp,'w') as dest:
            for k,v in source.attrs.items():dest.attrs[k]=v
            for k in source:source.copy(k,dest,name=k)
        temp.replace(path)
    finally:
        if temp.exists():temp.unlink()
    return {'path':str(path.relative_to(ROOT)),'before_bytes':before,'after_bytes':path.stat().st_size}

def main():
    compact=[compact_h5ad(p) for p in sorted(OUT.rglob('*.h5ad'))]
    obs,var,genes,x,raw=collapse_read()
    q=pd.read_csv(OUT/'snRNA/all_endothelial_QC_decisions.tsv.gz',sep='\t',index_col=0)
    ann=pd.read_csv(OUT/'snRNA/EC_annotations.tsv.gz',sep='\t',index_col=0)
    np.testing.assert_array_equal(q.index,obs.index)
    assert len(obs)==5573 and q.index.is_unique and ann.index.is_unique
    assert set(ann.index)==set(q.index[q.primary_EC_eligible])
    assert (~q.primary_EC_eligible).sum()+len(ann)==5573
    assert not q.loc[ann.index,['mixed_lineage_review','coherent_mural_group_review']].any().any()
    # Counts offered for download reconstruct the source for every retained AND reviewed nucleus.
    for stage in STAGES:
        for name,mask in [('EC_corrected_and_raw_counts.h5ad',q.primary_EC_eligible),('EC_review_corrected_and_raw_counts.h5ad',~q.primary_EC_eligible)]:
            expected=q.index[q.age.eq(stage)&mask]
            a=ad.read_h5ad(OUT/'snRNA'/stage/name)
            np.testing.assert_array_equal(a.obs_names,expected)
            np.testing.assert_array_equal(a.var_names,genes)
            ix=obs.index.get_indexer(expected)
            exact_sparse(a.X,x[ix]);exact_sparse(a.layers['counts_raw'],raw[ix])
            assert np.all(a.X.data>=0) and np.isfinite(a.X.data).all()
            assert np.allclose(a.layers['counts_raw'].data,np.round(a.layers['counts_raw'].data),atol=0,rtol=0)
    # The transfer's legacy shared-panel bundle contains no duplicate-symbol targets.
    with np.load(ROOT/'results/xenium_ventricular_mapping/reference_shared_panel_counts.npz') as b:
        ri=pd.Index(b['cell_ids'].astype(str)).get_indexer(obs.index);gi=pd.Index(genes).get_indexer(b['genes'].astype(str))
        np.testing.assert_allclose(x[:,gi].toarray(),b['corrected'][ri],rtol=2e-5,atol=1e-6)
        np.testing.assert_array_equal(raw[:,gi].toarray(),b['original'][ri])
    inherited=pd.read_csv(BASE/'results/snrna_qc/qc_summary.tsv',sep='\t')
    assert len(inherited)==8 and inherited.decontx_convergence.le(.001).all()
    inherited.to_csv(OUT/'snRNA/inherited_QC_summary.tsv',sep='\t',index=False)
    # Exact paper-panel coverage, including unavailable targets.
    coverage=[]
    with h5py.File(ROOT/'results/xenium_nuclear_anchored_segmentation/P0/nuclear_anchored_counts.h5ad') as h:
        spatial_genes=read_elem(h['var']).index.to_numpy()
    for family,panel in PANELS.items():
        for gene in panel:
            coverage.append({'family':family,'gene':gene,'RNA_measured':gene in genes,'Xenium_measured':gene in spatial_genes,'anchor':gene in ANCHORS[family]})
    save(pd.DataFrame(coverage),OUT/'subtype_panel_coverage.tsv')
    stages=[];allids=[]
    for stage in STAGES:
        folder=OUT/'Xenium'/stage
        a=pd.read_csv(folder/'EC_consensus_annotations.tsv.gz',sep='\t',index_col=0,float_precision='round_trip')
        atlas=pd.read_csv(ROOT/'results/whole_heart_spatial_atlas'/stage/'all_cell_annotations.tsv.gz',sep='\t',index_col=0,float_precision='round_trip')
        mask=(atlas.atlas_identity.eq('Endothelial')|atlas.whole_reference_lineage_proposal.eq('Endothelial')).to_numpy()
        np.testing.assert_array_equal(a.index,atlas.index[mask])
        np.testing.assert_array_equal(a[['x_centroid','y_centroid']].to_numpy(),atlas.loc[a.index,['x_centroid','y_centroid']].to_numpy())
        assert a.atlas_identity.equals(atlas.loc[a.index,'atlas_identity'])
        assert a.index.is_unique and atlas.loc[a.index,'pass_expression_QC'].all()
        cells=ad.read_h5ad(folder/'EC_candidate_count_views.h5ad')
        np.testing.assert_array_equal(cells.obs_names,a.index);np.testing.assert_array_equal(cells.var_names,spatial_genes)
        assert list(cells.obs)==list(a)
        assert cells.obs.confidence_tier.to_numpy().tolist()==a.confidence_tier.to_numpy().tolist()
        with h5py.File(ROOT/'results/xenium_nuclear_anchored_segmentation'/stage/'nuclear_anchored_counts.h5ad') as h:
            exact_sparse(cells.X,sparse_dataset(h['X'])[mask])
            for key in ['high_confidence','nuclear_only','depth_matched','nuclear_depth_matched']:
                exact_sparse(cells.layers[key],sparse_dataset(h['layers/'+key])[mask])
        with h5py.File(ROOT/'results/xenium_ventricular_mapping'/stage/f'{stage}_xenium_counts_and_candidate_mapping.h5ad') as h:
            exact_sparse(cells.layers['native'],sparse_dataset(h['X'])[mask])
        sums={k:np.asarray(v.sum(axis=1)).ravel() for k,v in [('new',cells.X)]+list(cells.layers.items())}
        np.testing.assert_array_equal(sums['nuclear_only'],sums['nuclear_depth_matched'])
        np.testing.assert_array_equal(sums['high_confidence'],sums['depth_matched'])
        assert (sums['nuclear_only']<=sums['new']).all() and (sums['high_confidence']<=sums['new']).all()
        cyto=cells.X-cells.layers['nuclear_only'];assert (cyto.data>=0).all()
        frames={}
        for v in ['native','nuclear_anchored','high_confidence','nuclear_only','nuclear_calibrated','depth_matched','nuclear_depth_matched','cytoplasmic_only']:
            f=pd.read_csv(folder/f'EC_{v}_assignments.tsv.gz',sep='\t',index_col=0)
            np.testing.assert_array_equal(f.index,a.index)
            counts=cells.X if v=='nuclear_anchored' else cyto if v=='cytoplasmic_only' else cells.layers['nuclear_only' if v=='nuclear_calibrated' else v]
            np.testing.assert_array_equal(np.asarray(counts.sum(axis=1)).ravel(),f.total_gene_molecules)
            # Verify markers in actual counts for every released individual-view call.
            for family in FAMILIES:
                m=f.fine_supported & f.mapped_family.eq(family)
                n,anchors,ok=marker_support(counts,spatial_genes,family)
                assert ok[m].all() and np.array_equal(n[m],f.loc[m,'marker_genes'])
                if m.any():
                    assert f.loc[m,'information_QC'].all() and f.loc[m,'broad_EC_competition_pass'].all()
                    assert f.loc[m,'both_library_agreement'].all() and f.loc[m,'raw_reference_agreement'].all()
                    assert not f.loc[m,'mixed_lineage_review'].any()
                    assert f.loc[m,'family_weight'].ge(.6).all() and f.loc[m,'family_margin'].ge(.15).all()
            frames[v]=f
        assert a.high_confidence_subtype.equals(a.confidence_tier.str.startswith(('A:','B:')))
        assert a.assigned_including_moderate.equals(a.confidence_tier.str.startswith(('A:','B:','C:')))
        for view in ['nuclear_anchored','nuclear_only','nuclear_calibrated']:
            m=a.high_confidence_subtype
            assert frames[view].loc[m,'fine_supported'].all()
            np.testing.assert_array_equal(a.loc[m,'fine_consensus_family'],frames[view].loc[m,'mapped_family'])
        for view in ['nuclear_only','nuclear_calibrated']:
            m=a.nuclear_supported_subtype
            assert frames[view].loc[m,'fine_supported'].all()
            np.testing.assert_array_equal(a.loc[m,'nuclear_supported_family'],frames[view].loc[m,'mapped_family'])
        assert not a.loc[a.multiple_projected_nuclear_components_review,['assigned_including_moderate','nuclear_supported_subtype']].any().any()
        assert not (a.nuclear_only_subtype_candidate & a.assigned_including_moderate).any()
        assert a.loc[~a.assigned_including_moderate,'fine_consensus_family'].eq(UNRESOLVED).all()
        assert not a.loc[a.high_confidence_subtype,'geometry_unresolved_review'].any()
        graph=pd.read_csv(folder/'EC_de_novo_clusters.tsv.gz',sep='\t',index_col=0)
        assert graph.index.isin(a.index).all()
        assert not a.loc[graph.index,'multiple_projected_nuclear_components_review'].any()
        for view in ['nuclear_anchored','nuclear_calibrated']:
            assert not frames[view].loc[graph.index,'mixed_lineage_review'].any()
        for view in ['nuclear_anchored','high_confidence','nuclear_only','nuclear_depth_matched']:
            assert frames[view].loc[graph.index,'information_QC'].all()
        # Every physical coordinate remains the native centroid; no smoothing or label quota.
        stages.append({'stage':stage,'candidate_profiles':len(a),'native_broad_EC':int(a.native_EC_broad.sum()),
            'high_confidence_subtypes':int(a.high_confidence_subtype.sum()),'extra_nuclear_candidates':int(a.nuclear_only_subtype_candidate.sum())})
        allids.extend(stage+'::'+i for i in a.index)
        print('VALIDATED_EC_STAGE',stage,flush=True)
    assert len(allids)==len(set(allids))==95316
    assert sum(s['native_broad_EC'] for s in stages)==69617
    for p in OUT.glob('*.png'):
        with Image.open(p) as im:im.verify()
    pdfpath=OUT/'endothelial_QC_and_subtypes_report.pdf'
    info=subprocess.check_output(['pdfinfo',str(pdfpath)],text=True)
    pages=int(re.search(r'^Pages:\s+(\d+)',info,re.MULTILINE).group(1));assert pages==27
    text=subprocess.check_output(['pdftotext',str(pdfpath),'-'],text=True)
    assert f'{len(ann):,}' in text and all(s in text for s in STAGES)
    assert len(pd.read_csv(OUT/'EC_interface_ROI_coverage.tsv',sep='\t'))==8
    assert all(p.stat().st_size<95*1024**2 for p in OUT.rglob('*') if p.is_file())
    validation={'status':'PASS','RNA_source_EC_nuclei':5573,'primary_RNA_nuclei':len(ann),'RNA_review_nuclei':int((~q.primary_EC_eligible).sum()),
        'RNA_raw_and_corrected_counts_exact_for_all_primary_and_review_nuclei':True,'inherited_eight_DecontX_fits_converged':True,
        'shared_panel_reference_counts_agree_with_full_gene_counts':True,'Xenium_stages':stages,
        'all_native_EC_profiles_accounted_for':69617,'all_EC_candidate_profiles_accounted_for':95316,
        'Xenium_native_new_core_nuclear_depth_counts_exactly_reproduce_existing_sources':True,
        'actual_measured_marker_support_validated_for_every_view_fine_label':True,'strict_and_extra_nuclear_support_rules_checked':True,
        'multiple_projected_nuclear_components_excluded_from_fine_release_not_deleted':True,
        'coordinates_and_parent_atlas_identities_unchanged':True,'PDF_pages':pages,'fixed_EC_interface_ROIs':8,
        'confidence_not_independent_biological_accuracy':True,'source_RNA_sha256':sha(RNA_INPUT)}
    validation['HDF5_free_space_compaction_then_exact_source_count_validation']=compact
    write_json(OUT/'validation.json',validation)
    paths=[p for p in OUT.rglob('*') if p.is_file() and p.name!='file_manifest.tsv']
    paths += [ROOT/'analysis'/s for s in ['endothelial_common.py','reassess_endothelial_RNA.py','reassess_endothelial_Xenium.py','consolidate_endothelial_QC.py','report_endothelial_QC.py','validate_endothelial_QC.py']]
    paths += [p for p in (ROOT/'analysis/endothelial_upstream_QC').glob('*') if p.is_file()]
    save(pd.DataFrame([{'path':str(p.relative_to(ROOT)),'bytes':p.stat().st_size,'sha256':sha(p)} for p in sorted(paths)]),OUT/'file_manifest.tsv')
    print(json.dumps(validation,indent=2),flush=True)

if __name__=='__main__':main()
