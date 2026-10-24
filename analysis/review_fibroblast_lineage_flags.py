"""Freeze residual lineage-mixing review flags before any spatial association fit."""
import anndata as ad
import numpy as np
import pandas as pd
from fibroblast_state_common import *

ann=pd.read_csv(OUT/'RNA_cell_annotations.tsv.gz',sep='\t',index_col=0)
scores=pd.read_csv(OUT/'RNA_cell_program_scores.tsv.gz',sep='\t',index_col=0)
summary=pd.read_csv(OUT/'RNA_cluster_summary.tsv',sep='\t')
means=ann[['age','fibroblast_cluster_id']].join(scores[['CM_carryover','mural_check','fibroblast_check']]).groupby('fibroblast_cluster_id',observed=True).mean(numeric_only=True)
rows=[];flags={}
for stage in STAGES:
    subset=ann.age.eq(stage)
    cm_median=float(scores.loc[subset,'CM_carryover'].median())
    mural_median=float(scores.loc[subset,'mural_check'].median())
    for cid in ann.loc[subset,'fibroblast_cluster_id'].unique():
        r=means.loc[cid]
        cmflag=r.CM_carryover>.75 and r.CM_carryover>2*cm_median
        muralflag=r.mural_check>.3 and r.mural_check>3*mural_median
        reason='CM-enriched residual program' if cmflag else 'mural-enriched residual program' if muralflag else 'None'
        flags[cid]=bool(cmflag or muralflag)
        rows.append({'stage':stage,'cluster_id':cid,'CM_mean_score':r.CM_carryover,'stage_CM_median':cm_median,
            'mural_mean_score':r.mural_check,'stage_mural_median':mural_median,
            'fibroblast_mean_score':r.fibroblast_check,'lineage_mixing_review_flag':flags[cid],'reason':reason})
review=pd.DataFrame(rows)
save_tsv(review,OUT/'RNA_residual_lineage_review.tsv')
if 'family_before_lineage_review' not in ann:ann['family_before_lineage_review']=ann.fibroblast_family
ann['lineage_mixing_review_flag']=ann.fibroblast_cluster_id.map(flags)
ann['fine_reference_eligible']=~ann.lineage_mixing_review_flag
ann.loc[ann.lineage_mixing_review_flag,'fibroblast_family']=UNRESOLVED
ann.loc[ann.lineage_mixing_review_flag,'fibroblast_support']='lineage-mixing review; withheld from fine reference'
ann.to_csv(OUT/'RNA_cell_annotations.tsv.gz',sep='\t',na_rep='NA')
if 'family_before_lineage_review' not in summary:summary['family_before_lineage_review']=summary.family
summary['lineage_mixing_review_flag']=summary.cluster_id.map(flags)
summary.loc[summary.lineage_mixing_review_flag,'family']=UNRESOLVED
summary.loc[summary.lineage_mixing_review_flag,'support']='lineage-mixing review; withheld from fine reference'
save_tsv(summary,OUT/'RNA_cluster_summary.tsv')
for stage in STAGES:
    part=ann[ann.age.eq(stage)]
    part.to_csv(OUT/stage/'RNA_cell_annotations.tsv.gz',sep='\t',na_rep='NA')
    p=OUT/stage/'RNA_embeddings_metadata_only.h5ad'
    obj=ad.read_h5ad(p);obj.obs=part.loc[obj.obs_names].copy();obj.write_h5ad(p,compression='gzip')
meta=pd.read_csv(OUT/'RNA_pseudobulk_metadata.tsv',sep='\t',index_col=0)
mapper=ann.drop_duplicates('fibroblast_cluster_id').set_index('fibroblast_cluster_id')
meta['fibroblast_family']=meta.fibroblast_cluster_id.map(mapper.fibroblast_family)
meta['fine_reference_eligible']=~meta.fibroblast_cluster_id.map(flags)
meta.to_csv(OUT/'RNA_pseudobulk_metadata.tsv',sep='\t')
combined=ann[['sample','age','fibroblast_family','fine_reference_eligible','total_counts','decontx_contamination']].join(scores)
group=combined.groupby(['sample','age','fibroblast_family','fine_reference_eligible'],observed=True)
unit=group[list(scores)].mean().join(group.size().rename('nuclei')).reset_index()
save_tsv(unit,OUT/'RNA_matrix_programs_by_library_family.tsv')
write_json(OUT/'RNA_lineage_review_methods.json',{'frozen_before_spatial_models':True,
    'CM_flag':'Cluster mean CM-specific panel score>0.75 and>2x same-stage median across nuclei',
    'mural_flag':'Cluster mean mural-panel score>0.3 and>3x same-stage median',
    'purpose':'Conservative residual lineage-mixing review, not validated doublet/contamination rate',
    'nuclei_flagged':int(ann.lineage_mixing_review_flag.sum()),'flagged_clusters':review.loc[review.lineage_mixing_review_flag,'cluster_id'].tolist(),
    'handling':'All source nuclei/counts retained, label set unresolved, excluded from fine reference and primary matrix pseudobulk comparisons; original family proposal retained'})
print(review[review.lineage_mixing_review_flag].to_string(index=False),flush=True)
