"""Test identity and Ankrd1 sensitivity to measured molecule assignments."""
import argparse
import json

import numpy as np
import pandas as pd

from design_xenium_boundary_qc import OUT,ATLAS,FIT,STAGES
from process_xenium_boundary_molecules import VARIANTS
from boundary_reference_projection import FrozenProjection
from reassess_whole_heart_lineages import PANELS,WITHHOLD,digest


def pair_effects(stage,variant,scope,ids,counts,genes,obs,classification,roi_by_cell=None):
    pairs=pd.read_csv(FIT/stage/'default_contrast_matched_pairs.tsv.gz',sep='\t')
    pairs=pairs[pairs.definition.eq('independent_strict')&pairs.cohort.eq('whole_reference_strict')].copy()
    index=pd.Index(ids);ni=index.get_indexer(pairs.near_cell_id);fi=index.get_indexer(pairs.far_cell_id)
    ok=(ni>=0)&(fi>=0)
    if roi_by_cell is not None:
        both=ok.copy();ok[both]=np.array([roi_by_cell[pairs.near_cell_id.iloc[i]]==roi_by_cell[pairs.far_cell_id.iloc[i]] for i in np.flatnonzero(both)])
    pairs=pairs.loc[ok].copy();ni=ni[ok];fi=fi[ok]
    ank=counts[:,list(genes).index('Ankrd1')].astype(float)
    denominator=counts.sum(1)-counts[:,np.isin(genes,WITHHOLD)].sum(1)
    expression=1e4*ank/np.maximum(denominator,1);log=np.log1p(expression)
    nuclear=obs.loc[ids,'nucleus_area'].isna().to_numpy();missing=nuclear[ni].astype(int)+nuclear[fi].astype(int)
    passing=classification.information_QC_pass.to_numpy()
    cm=classification.new_identity.eq('Cardiomyocyte').to_numpy()
    rows=[]
    for cohort,keep in [('all_frozen_pairs',np.ones(len(pairs),bool)),('both_information_QC',passing[ni]&passing[fi]),
            ('both_CM_identity_retained',passing[ni]&passing[fi]&cm[ni]&cm[fi])]:
        for stratum,selected in [('all',keep),('both_nuclear_area_measured',keep&(missing==0)),
                ('one_missing_nuclear_area',keep&(missing==1)),('both_missing_nuclear_area',keep&(missing==2))]:
            near,far=ni[selected],fi[selected];segments=pairs.segment_id.to_numpy()[selected]
            if len(near):
                segment=pd.Series(log[near]-log[far]).groupby(segments).mean()
                ratio=expression[near].mean()/max(expression[far].mean(),1e-12)
                mean=float(segment.mean());positive=float(segment.gt(0).mean())
            else:ratio=mean=positive=np.nan
            rows.append({'stage':stage,'variant':variant,'scope':scope,'pair_filter':cohort,'nuclear_area_stratum':stratum,
                'pairs':len(near),'segments':len(set(segments)),'near_far_mean_Ankrd1_CP10k_ratio':ratio,
                'equal_segment_mean_log_difference':mean,'fraction_segments_positive':positive,
                'near_mean_gene_transcripts':float(counts[near].sum(1).mean()) if len(near) else np.nan,
                'far_mean_gene_transcripts':float(counts[far].sum(1).mean()) if len(far) else np.nan})
    return rows


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--stage',choices=STAGES,required=True);args=parser.parse_args();stage=args.stage
    folder=OUT/stage
    if (folder/'sensitivity_completion.json').exists():print('SKIP_SENSITIVITY',stage,flush=True);return
    assert (folder/'ROI_reassignment_completion.json').exists()
    obs=pd.read_csv(ATLAS/stage/'all_cell_annotations.tsv.gz',sep='\t',index_col=0,float_precision='round_trip')
    target=pd.read_csv(folder/'boundary_target_cells.tsv.gz',sep='\t',index_col=0,float_precision='round_trip')
    with np.load(folder/'geometry_sensitivity_counts.npz') as bundle:
        ids,genes=bundle['cell_ids'].astype(str),bundle['genes'].astype(str)
        matrices={k:bundle[k] for k in VARIANTS};regions=bundle['region_counts'];region_names=bundle['region_names'].astype(str)
    np.testing.assert_array_equal(ids,target.index)
    # Marker localization is measured directly, independently of any variant fit.
    counter=np.zeros(regions.shape[:2],np.int64);own=np.zeros_like(counter)
    for label in target.audit_lineage.unique():
        selected=target.audit_lineage.eq(label).to_numpy();panel=PANELS.get(label,[])
        foreign=set().union(*(set(PANELS[k]) for k in ['Endothelial','Fibroblast','Pericyte','Smooth muscle','Immune'])) if label=='Cardiomyocyte' else set(PANELS['Cardiomyocyte'])
        counter[:,selected]=regions[:,selected,:][:,:,np.isin(genes,list(foreign))].sum(2)
        own[:,selected]=regions[:,selected,:][:,:,np.isin(genes,panel)].sum(2)
    totals=regions.sum(2);marker_rows=[]
    for j,name in enumerate(region_names):
        part=target[['atlas_identity','audit_lineage','nuclear_area_missing','ROI_inner_evaluation_cell']].copy()
        part['region']=name;part['gene_transcripts']=totals[j];part['counter_lineage_marker_transcripts']=counter[j]
        part['own_lineage_marker_transcripts']=own[j]
        part['counter_marker_fraction']=np.divide(counter[j],totals[j],out=np.full(len(target),np.nan),where=totals[j]>0)
        part['own_marker_fraction']=np.divide(own[j],totals[j],out=np.full(len(target),np.nan),where=totals[j]>0)
        marker_rows.append(part.reset_index())
    markers=pd.concat(marker_rows,ignore_index=True);markers.to_csv(folder/'cell_shell_marker_profiles.tsv.gz',sep='\t',index=False,na_rep='NA')
    group=markers.groupby(['atlas_identity','audit_lineage','region'],dropna=False).agg(cells_with_molecules=('gene_transcripts',lambda x:int(x.gt(0).sum())),
        gene_transcripts=('gene_transcripts','sum'),counter_marker_transcripts=('counter_lineage_marker_transcripts','sum'),
        own_marker_transcripts=('own_lineage_marker_transcripts','sum')).reset_index()
    group['counter_marker_fraction']=group.counter_marker_transcripts/np.maximum(group.gene_transcripts,1)
    group['own_marker_fraction']=group.own_marker_transcripts/np.maximum(group.gene_transcripts,1)
    group['stage']=stage;group.to_csv(folder/'boundary_marker_localization_summary.tsv',sep='\t',index=False,na_rep='NA')
    # Whole-run quality uses the correct gene-only denominator, matching vendor metrics.
    gene_totals=pd.read_csv(folder/'molecule_gene_totals.tsv',sep='\t')
    run=json.loads((folder/'molecule_run_metrics.json').read_text())
    run['observed_gene_only_Q20_fraction']=float(gene_totals.QV20.sum()/gene_totals.all_QV.sum())
    assert abs(run['observed_gene_only_Q20_fraction']-run['vendor_Q20_fraction'])<1e-10
    assert abs(run['observed_Q20_gene_assigned_fraction']-run['vendor_gene_assigned_fraction'])<1e-10
    (folder/'verified_gene_QV_metrics.json').write_text(json.dumps(run,indent=2)+'\n')
    model=FrozenProjection(stage,obs,genes)
    np.savez_compressed(folder/'frozen_query_calibration.npz',query_mean=model.query_mean,query_sd=model.query_sd,reference_ids=model.refids)
    effects=[];retention=[];assignment_paths=[];all_transitions=[]

    def evaluate(variant,scope,cell_ids,counts,baseline_counts,roi_by_cell=None):
        dest=folder/(variant+'_projection_assignments.tsv.gz')
        if dest.exists():
            classified=pd.read_csv(dest,sep='\t',index_col=0,float_precision='round_trip')
            np.testing.assert_array_equal(classified.index,cell_ids)
        else:
            classified=model.classify(cell_ids,counts);classified.to_csv(dest,sep='\t',na_rep='NA')
        assignment_paths.append(dest.name)
        if variant=='native_Q20':
            expected=obs.loc[cell_ids,'atlas_identity'];difference=classified.new_identity.ne(expected)
            assert not difference.any(),classified.loc[difference].join(expected.rename('expected')).to_dict('index')
            (folder/'native_projection_validation.json').write_text(json.dumps({'stage':stage,'cells':len(cell_ids),
                'exact_existing_atlas_identity_reproduction':True,'mismatches':0},indent=2)+'\n')
        effects.extend(pair_effects(stage,variant,scope,cell_ids,counts,genes,obs,classified,roi_by_cell))
        baseline=obs.loc[cell_ids,'atlas_identity'];nuc=obs.loc[cell_ids,'nucleus_area'].isna().to_numpy()
        for label in sorted(baseline.unique()):
            for nuclear_group,flag in [('all',np.ones(len(cell_ids),bool)),('measured',~nuc),('missing',nuc)]:
                select=baseline.eq(label).to_numpy()&flag;n=int(select.sum())
                if not n:continue
                before=baseline_counts[select].sum(1);after=counts[select].sum(1)
                retention.append({'stage':stage,'variant':variant,'scope':scope,'baseline_identity':label,'nuclear_area':nuclear_group,
                    'profiles':n,'information_QC_pass':int(classified.information_QC_pass.to_numpy()[select].sum()),
                    'identity_unchanged':int(np.equal(classified.new_identity.to_numpy()[select],baseline.to_numpy()[select]).sum()),
                    'mean_transcript_retained_ratio':float(np.mean(after/np.maximum(before,1))),
                    'median_transcript_retained_ratio':float(np.median(after/np.maximum(before,1)))})
        transition=pd.crosstab(baseline,classified.new_identity)
        for old,row in transition.iterrows():
            for new,n in row.items():
                if n:all_transitions.append({'stage':stage,'variant':variant,'scope':scope,'original_identity':old,'new_identity':new,'profiles':int(n)})
        print('CLASSIFIED',stage,variant,scope,len(cell_ids),flush=True)

    for variant in VARIANTS:evaluate(variant,'geometry_targets_and_all_frozen_pairs',ids,matrices[variant],matrices['native_Q20'])
    # Unique ownership prevents cells in overlapping image crops being counted twice.
    rois=pd.read_csv(OUT/'selected_ROIs.tsv',sep='\t');rois=rois[rois.stage.eq(stage)]
    seen=set();morph_ids=[];roi_by_cell={};collected={v:[] for v in ['native_raster_Q20','centroid_Voronoi','membrane_watershed']}
    baseline=[]
    for roi in rois.itertuples():
        with np.load(folder/(roi.roi_id+'_ROI_count_variants.npz')) as b:
            cells=b['cell_ids'].astype(str);mask=b['comparison_eligible']&~np.isin(cells,list(seen));picked=cells[mask]
            ix=pd.Index(ids).get_indexer(picked);assert (ix>=0).all()
            baseline.append(matrices['native_Q20'][ix]);morph_ids.extend(picked);seen.update(picked)
            roi_by_cell.update({c:roi.roi_id for c in picked})
            for variant in collected:collected[variant].append(b[variant][mask])
    morph_ids=np.array(morph_ids,dtype=str);baseline=np.concatenate(baseline)
    morphology={v:np.concatenate(parts) for v,parts in collected.items()}
    np.savez_compressed(folder/'morphology_comparison_counts.npz',cell_ids=morph_ids,genes=genes,
        native_Q20=baseline,roi_ids=np.array([roi_by_cell[c] for c in morph_ids]),**morphology)
    # Fair source baseline on this smaller image-covered cohort uses original molecules.
    native=pd.read_csv(folder/'native_Q20_projection_assignments.tsv.gz',sep='\t',index_col=0,float_precision='round_trip').loc[morph_ids]
    effects.extend(pair_effects(stage,'native_Q20','image_ROI_subset',morph_ids,baseline,genes,obs,native,roi_by_cell))
    for variant,counts in morphology.items():evaluate(variant,'image_ROI_subset',morph_ids,counts,baseline,roi_by_cell)
    pd.DataFrame(effects).to_csv(folder/'Ankrd1_boundary_sensitivity.tsv',sep='\t',index=False,na_rep='NA')
    pd.DataFrame(retention).to_csv(folder/'identity_and_information_retention.tsv',sep='\t',index=False,na_rep='NA')
    pd.DataFrame(all_transitions).to_csv(folder/'boundary_identity_transitions.tsv',sep='\t',index=False,na_rep='NA')
    (folder/'sensitivity_completion.json').write_text(json.dumps({'stage':stage,'geometry_targets':len(ids),'image_comparison_cells':len(morph_ids),
        'assignment_tables':assignment_paths,'frozen_reference_and_query_calibration':True,'nucleus_required':False,
        'ROIs_deduplicated_for_cell_summary':True,'all_source_baseline_labels_reproduced':True,
        'molecule_completion_sha256':digest(folder/'molecule_completion.json'),
        'morphology_completion_sha256':digest(folder/'morphology_completion.json'),
        'ROI_reassignment_completion_sha256':digest(folder/'ROI_reassignment_completion.json')},indent=2)+'\n')
    print('BOUNDARY_SENSITIVITY_COMPLETE',stage,flush=True)


if __name__=='__main__':main()
