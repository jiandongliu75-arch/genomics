"""Audit fibroblast attrition and add nuclear-optional spatial sensitivities.

Frozen atlas labels and previous releases are unchanged. A multigene candidate
layer additionally removes the compulsory single-anchor gate, while preserving
reference, cell/core, marker-count and positive nuclear-counterevidence checks.
It is a sensitivity layer, not a newly validated cell identity atlas.
"""
from pathlib import Path
import os,json,hashlib,argparse
os.environ.setdefault('OPENBLAS_NUM_THREADS','2')
import numpy as np
import pandas as pd
from scipy import sparse
import cardiac_subtype_niches as old

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/cardiac_biological_synthesis'
STAGES=old.STAGES;FB=old.FB
INPUTS=set()
def read(path,**kw):
    path=Path(path);INPUTS.add(str(path.relative_to(ROOT)))
    return pd.read_csv(path,sep='\t',**kw)
def save(d,name):
    OUT.mkdir(exist_ok=True,parents=True)
    d.to_csv(OUT/name,sep='\t',index=False,na_rep='NA',float_format='%.9g')
def gates(d,anchor=True):
    ok=(d.family_proposal.isin(FB)&d.family_weight.ge(.6)&d.family_margin.ge(.15)
        &d.both_library_agreement&d.raw_reference_agreement
        &d.measured_state_markers.ge(2)&d.state_information_QC
        &~d.mixed_CM_FB_marker_review)
    return ok&d.measured_state_anchors.ge(1) if anchor else ok

def annotations(stage):
    folder=ROOT/'results/xenium_nuclear_anchored_segmentation'/stage
    v={name:read(folder/f'fibroblast_{name}_assignments.tsv.gz',index_col=0)
       for name in ['native','nuclear_anchored','high_confidence','nuclear_only','nuclear_modality_calibrated']}
    a,c,n,k=[v[t] for t in ['nuclear_anchored','high_confidence','nuclear_only','nuclear_modality_calibrated']]
    atlas=read(ROOT/'results/cardiac_spatial_pipeline'/stage/'atlas_annotations.tsv.gz',index_col=0)
    f=atlas.reindex(a.index);assert all(x.index.equals(a.index) for x in v.values())
    eligible=f.atlas_identity.eq('Fibroblast')&f.program_source_eligible
    same=a.family_proposal.eq(c.family_proposal)
    # A supported alternative nuclear identity is counterevidence. Missing or
    # insufficient nuclear support never vetoes a concordant cell/core call.
    veto=(n.supported_family_assignment&n.family_proposal.ne(a.family_proposal))|(
        k.supported_family_assignment&k.family_proposal.ne(a.family_proposal))
    whole=eligible&gates(a)
    boundary=whole&gates(c)&same&~veto
    program=eligible&gates(a,False)&gates(c,False)&same&~veto
    legacy=eligible&f.fine_identity.isin(FB)
    out=f[['stage','x_centroid','y_centroid','atlas_identity','author_atlas_cell_state','fine_identity','program_source_eligible']].copy()
    out['legacy_strict_family']=f.fine_identity.where(legacy,'Unresolved')
    for key,mask in [('whole_cell_supported',whole),('boundary_supported_nuclear_optional',boundary),('multigene_boundary_candidate',program)]:
        out[key]=a.family_proposal.where(mask,'Unresolved')
    out['supported_nuclear_counterevidence']=veto
    for key in ['family_proposal','family_weight','family_margin','both_library_agreement','raw_reference_agreement',
                'measured_state_markers','measured_state_anchors','state_information_QC','mixed_CM_FB_marker_review']:
        out['whole_cell__'+key]=a[key]
    out['core_family_proposal']=c.family_proposal
    out['core_markers']=c.measured_state_markers
    out['core_anchors']=c.measured_state_anchors
    out['nuclear_supported_family']=n.family_proposal.where(n.supported_family_assignment,'Unresolved')
    out.index.name='cell_id'
    return atlas,out,v,eligible

def prepare():
    coverage=[];attrition=[];confusions=[];postn=[]
    for stage in STAGES:
        print('FB_SENSITIVITY_PREPARE',stage,flush=True)
        atlas,d,v,eligible=annotations(stage);a=v['nuclear_anchored']
        for family in FB:
            pub=family.split('-')[0]+'.Fib'
            row={'stage':stage,'family':family,'eligible_FB_profiles':int(eligible.sum()),
                 'public_comparator':int((eligible&d.author_atlas_cell_state.eq(pub)).sum())}
            for col in ['legacy_strict_family','whole_cell_supported','boundary_supported_nuclear_optional','multigene_boundary_candidate']:
                row[col]=int(d[col].eq(family).sum())
            coverage.append(row)
        for pub,g in d.loc[eligible].groupby('author_atlas_cell_state'):
            z=a.loc[g.index]
            for proposal,count in z.family_proposal.value_counts().items():
                confusions.append({'stage':stage,'public_label':pub,'whole_cell_proposal':proposal,'profiles':int(count)})
            expected=pub.split('.')[0]+'-associated'
            if expected not in FB:continue
            checks={'matching_whole_cell_proposal':z.family_proposal.eq(expected),'weight_ge_0.6':z.family_weight.ge(.6),
                    'margin_ge_0.15':z.family_margin.ge(.15),'both_reference_libraries':z.both_library_agreement,
                    'raw_reference_agreement':z.raw_reference_agreement,'two_measured_markers':z.measured_state_markers.ge(2),
                    'one_anchor':z.measured_state_anchors.ge(1),'capture_information_QC':z.state_information_QC}
            for gate,passed in checks.items():
                attrition.append({'stage':stage,'public_label':pub,'gate':gate,'profiles':len(z),'pass_gate_independently':int(passed.sum()),
                    'interpretation':'Overlapping gates, not sequential losses; not classification accuracy'})
        save(d.reset_index(),stage+'_FB_assignment_evidence.tsv.gz')
        # Direct measured Postn evidence; no reference-label requirement.
        import anndata as ad
        path=ROOT/'results/xenium_ventricular_mapping'/stage/(stage+'_xenium_counts_and_candidate_mapping.h5ad')
        INPUTS.add(str(path.relative_to(ROOT)));b=ad.read_h5ad(path,backed='r')
        ids=d.index[eligible];ix=b.obs_names.get_indexer(ids);assert (ix>=0).all()
        x=b.X[ix].tocsr();gi=np.flatnonzero(b.var_names=='Postn');assert len(gi)==1
        y=x[:,gi[0]].toarray().ravel();den=np.asarray(x.sum(1)).ravel()
        postn.append({'stage':stage,'eligible_FB_profiles':len(ids),'Postn_detected_profiles':int((y>0).sum()),
            'Postn_detection_fraction':float((y>0).mean()),'Postn_mean_log1p_per10000_panel_molecules':float(np.log1p(y*1e4/np.maximum(den,1)).mean()),
            'Postn_total_molecules':float(y.sum()),'panel_capture_total':float(den.sum()),
            'scope':'Direct native Xenium detection within eligible broad FBs, not subtype abundance or full-transcriptome normalization'})
        b.file.close()
    save(pd.DataFrame(coverage),'FB_assignment_coverage_comparison.tsv')
    save(pd.DataFrame(attrition),'FB_assignment_gate_audit.tsv')
    save(pd.DataFrame(confusions),'FB_public_vs_whole_cell_proposals.tsv')
    save(pd.DataFrame(postn),'Postn_direct_Xenium_evidence.tsv')
    markers=read(ROOT/'results/cardiac_developmental_story/fibroblast_library_marker_profiles.tsv')
    save(markers[markers.level.eq('fibroblast_family')&markers.gene.eq('Postn')],'Postn_direct_snRNA_evidence.tsv')
    (OUT/'assignment_rules.json').write_text(json.dumps({
        'legacy_atlas_unchanged':True,'nuclear_optional':True,
        'whole_cell':'Inherited whole-cell support gates, eligible broad FB and segmentation/source QC',
        'boundary_supported':'Whole-cell and boundary-core agree and pass inherited gates; supported alternative nuclear call vetoes; nuclear absence never vetoes',
        'multigene_candidate':'Same as boundary-supported, except no mandatory individual anchor in either cell/core; >=2 supplied markers, weight>=0.6, margin>=0.15, both RNA libraries, raw-reference agreement and information QC retained in both views',
        'nuclear_evidence_is_independent_validation':False,'reference_weights_are_calibrated_probabilities':False,
        'candidate_layer_is_validated_identity':False,'spatial_smoothing_or_label_quotas':False,
        'no_cross_stage_Postn_label_restoration':True,'boundary_core_and_whole_cell_share_transcripts':True},indent=2)+'\n')
    (OUT/'audit_input_paths.json').write_text(json.dumps(sorted(INPUTS),indent=2)+'\n')

def models():
    allrows=[];contrasts=[]
    cov=read(OUT/'FB_assignment_coverage_comparison.tsv')
    original_coverage=read(ROOT/'results/cardiac_subtype_niches/source_subtype_coverage.tsv')
    for stage in STAGES:
        atlas=read(ROOT/'results/cardiac_spatial_pipeline'/stage/'atlas_annotations.tsv.gz',index_col=0)
        fb=read(OUT/(stage+'_FB_assignment_evidence.tsv.gz'),index_col='cell_id')
        for geometry in ['proximity50','proximity30','contact']:
            path=ROOT/'results/cardiac_spatial_pipeline'/stage/(geometry+'_graph.npz');INPUTS.add(str(path.relative_to(ROOT)))
            graph=sparse.load_npz(path).tocsr()
            d=read(ROOT/'results/cardiac_subtype_niches'/(stage+'_'+geometry+'_CM_subtype_neighborhoods.tsv.gz'),index_col='cell_id')
            d=d[d.analysis_eligible].copy();g=graph[atlas.index.get_indexer(d.index)]
            for view in ['boundary_supported_nuclear_optional','multigene_boundary_candidate']:
                if geometry!='proximity50' and view!='multigene_boundary_candidate':continue
                print('FB_SENSITIVITY_MODEL',stage,geometry,view,flush=True)
                frame=d.copy();ex=[]
                for family in FB:
                    mask=pd.Series(False,index=atlas.index);mask.loc[fb.index]=fb[view].eq(family).to_numpy()
                    key='FB__'+family
                    frame[key]=np.asarray(g@mask.to_numpy(float)).ravel()/np.maximum(frame.eligible_FB_neighbors.to_numpy(),1)
                    n=int(cov.loc[cov.stage.eq(stage)&cov.family.eq(family),view].iloc[0])
                    if n>=30 and frame[key].gt(0).sum()>=50:ex.append(key)
                for key in original_coverage.loc[original_coverage.stage.eq(stage)&original_coverage.subtype.str.startswith('EC__')&original_coverage.primary_source_profiles.ge(30),'subtype']:
                    if frame[key].gt(0).sum()>=50:ex.append(key)
                result,delta=old.fit_states(frame,ex,geometry,'primary','within_region')
                for r in result:r['source_view']=view
                for r in delta:r['source_view']=view
                allrows.extend(result);contrasts.extend(delta)
                if geometry=='proximity50' and view=='multigene_boundary_candidate':
                    keep=['stage','CM_state','x_centroid','y_centroid','maturation_group','core_Z__maturation21','core_Z__metabolic','core_Z__sarcomere',
                          'eligible_FB_neighbors','eligible_EC_neighbors']+['FB__'+family for family in FB]
                    save(frame[keep].reset_index(),stage+'_candidate_FB_CM_neighborhoods.tsv.gz')
    out=pd.DataFrame(allrows);out['within_section_BH_q']=np.nan
    for _,g in out.groupby(['stage','geometry','source_view','adjustment']):
        p=g.within_section_approx_p.to_numpy();order=np.argsort(p);q=np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]
        out.loc[g.index[order],'within_section_BH_q']=np.minimum(q,1)
    save(out,'nuclear_optional_FB_maturation_coefficients.tsv')
    save(pd.DataFrame(contrasts),'nuclear_optional_FB_CM_identity_slope_differences.tsv')
    inputs=set(json.loads((OUT/'audit_input_paths.json').read_text()))|INPUTS
    (OUT/'audit_input_paths.json').write_text(json.dumps(sorted(inputs),indent=2)+'\n')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--step',choices=['prepare','models','all'],default='all');a=p.parse_args()
    if a.step in ['prepare','all']:prepare()
    if a.step in ['models','all']:models()
