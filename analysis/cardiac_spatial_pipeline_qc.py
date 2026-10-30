"""All-lineage ownership audit and a separately identified conservative sensitivity."""
from cardiac_spatial_pipeline import *
from cardiac_spatial_models import design,exposures,joint_fit,OUTCOMES

def audit():
    summaries=[]
    for stage in STAGES:
        print('ALL_LINEAGE_QC',stage,flush=True)
        f=frame(stage)
        b=ad.read_h5ad(ROOT/f'results/xenium_nuclear_anchored_segmentation/{stage}/nuclear_anchored_counts.h5ad')
        genes=b.var_names.to_numpy(dtype=str);core=b.layers['high_confidence'].tocsr();nuc=b.layers['nuclear_only'].tocsr()
        q=f[['stage_cell_id','atlas_identity','x_centroid','y_centroid']].copy()
        for view,x in [('core',core),('nuclear',nuc)]:
            total=np.asarray(x.sum(1)).ravel();strong={}
            for l in LINEAGES:
                ix=np.isin(genes,MARKERS[l]);mols=np.asarray(x[:,ix].sum(1)).ravel();det=np.asarray((x[:,ix]>0).sum(1)).ravel()
                q[view+'__'+l+'__genes']=det;q[view+'__'+l+'__fraction']=mols/np.maximum(total,1)
                strong[l]=(det>=3)&(mols>=10)&(mols/np.maximum(total,1)>=.10)
            q[view+'_strong_lineage_panels']=np.sum(list(strong.values()),axis=0)
        own=np.array([q['core__'+l+'__genes'].iloc[i] if l in LINEAGES else 0 for i,l in enumerate(f.atlas_identity)])
        need=np.where(f.atlas_identity.eq('Epicardial'),1,2)
        q['weak_core_own_lineage']=f.atlas_resolved&(own<need)
        # Related mural/stromal programs can overlap biologically; this is review,
        # not an automatic doublet removal or a deconvolution purity estimate.
        q['multiple_strong_core_lineages']=f.atlas_resolved&q.core_strong_lineage_panels.ge(2)
        q['multiple_strong_nuclear_lineages']=f.atlas_resolved&q.nuclear_strong_lineage_panels.ge(2)
        q['conservative_QC_review']=q.weak_core_own_lineage|q.multiple_strong_core_lineages|q.multiple_strong_nuclear_lineages|f.geometry_unseeded
        q['nucleus_missing']=f.nucleus_area.isna()
        q['primary_atlas_identity_preserved']=True
        save(q,f'{stage}/all_lineage_count_view_QC.tsv.gz',index=True)
        for l,z in q.groupby('atlas_identity',observed=True):
            summaries.append({'stage':stage,'identity':l,'profiles':len(z),'weak_core_own_lineage':int(z.weak_core_own_lineage.sum()),
                'multiple_strong_core_lineages':int(z.multiple_strong_core_lineages.sum()),'multiple_strong_nuclear_lineages':int(z.multiple_strong_nuclear_lineages.sum()),
                'conservative_QC_review':int(z.conservative_QC_review.sum()),'missing_nucleus_not_automatically_removed':int(z.nucleus_missing.sum())})
        del b,core,nuc;gc.collect()
    save(pd.DataFrame(summaries),'all_lineage_QC_summary.tsv')

def masked_frame(stage):
    f=frame(stage).copy();q=read(f"{CFG['output']}/{stage}/all_lineage_count_view_QC.tsv.gz")
    bad=q.conservative_QC_review.to_numpy(bool)
    # This is a sensitivity view; the released source atlas_identity stays intact.
    f.loc[bad,'atlas_identity']='Ambiguous/review';f.loc[bad,'atlas_resolved']=False
    f.loc[bad,'program_source_eligible']=False;f.loc[bad,PROGS]=np.nan
    return f,q

def sensitivities():
    from cardiac_spatial_pipeline import neighborhood
    rows=[];coverage=[]
    for stage in STAGES:
        print('QC_SENSITIVITY',stage,flush=True)
        f,q=masked_frame(stage)
        # Inject only the separately named diagnostic frame into neighborhood()
        # while preserving all source files and baseline results.
        import cardiac_spatial_pipeline as coremodule
        original=coremodule.frame
        try:
            coremodule.frame=lambda s:f if s==stage else original(s)
            d=neighborhood(stage,'proximity50')
        finally:coremodule.frame=original
        baseline=read(f"{CFG['output']}/{stage}/proximity50_CM_neighborhoods.tsv.gz")
        ids=baseline.index[baseline.geometry_eligible]
        good=~q.loc[ids,'conservative_QC_review']
        good&=d.loc[ids,'resolved_neighbors'].ge(CFG['model_minimum_neighbors'])&d.loc[ids,'Endothelial_eligible_neighbors'].ge(1)&d.loc[ids,'Fibroblast_eligible_neighbors'].ge(1)
        ids=ids[good];c=frame(stage).loc[ids];data=d.loc[ids]
        coverage.append({'stage':stage,'baseline_CMs':int(baseline.geometry_eligible.sum()),'conservative_masked_CMs':len(ids)})
        save(data,f'{stage}/QC_masked_proximity50_CM_neighborhoods.tsv.gz',index=True)
        if len(ids)<200:continue
        ee=exposures(data);ee=ee[[key for key in ee if ee[key].std()>1e-8 and (ee[key].abs()>1e-8).sum()>=50]]
        # Identical focal CMs isolate the count-view/source-mask sensitivity.
        for source_name,env in [('masked_sources',data),('baseline_sources_same_CMs',baseline.loc[ids])]:
            e=exposures(env)[list(ee)]
            for level in ['technical','anatomy','CM_state']:
                for r in joint_fit(c[OUTCOMES],design(c,env,level),e,c[['x_centroid','y_centroid']].to_numpy(),CFG['bootstrap_sensitivity']):
                    rows.append({'stage':stage,'source_view':source_name,'adjustment':level,**r})
    save(pd.DataFrame(rows),'conservative_QC_association_sensitivity.tsv.gz');save(pd.DataFrame(coverage),'conservative_QC_model_coverage.tsv')
    finalize_priorities()

def finalize_priorities():
    candidates=pd.read_csv(OUT/'candidate_relationships.tsv',sep='\t')
    qc=pd.read_csv(OUT/'conservative_QC_association_sensitivity.tsv.gz',sep='\t')
    candidates['pre_QC_priority_rule_passed']=candidates.get('pre_QC_priority_rule_passed',candidates.priority_within_section_candidate)
    candidates['QC_masked_coefficient']=np.nan;candidates['QC_masked_same_sign']=False;candidates['QC_masked_interval_same_direction']=False
    for i,r in candidates.iterrows():
        z=qc.loc[qc.stage.eq(r.stage)&qc.exposure.eq(r.exposure)&qc.adjustment.eq('anatomy')&qc.outcome.eq('CM_maturation_core')&qc.source_view.eq('masked_sources')]
        if len(z)==1:
            v=z.iloc[0];candidates.loc[i,'QC_masked_coefficient']=v.standardized_coefficient
            candidates.loc[i,'QC_masked_same_sign']=np.sign(v.standardized_coefficient)==np.sign(r.standardized_coefficient)
            candidates.loc[i,'QC_masked_interval_same_direction']=(v.block_CI_low>0 and r.standardized_coefficient>0) or (v.block_CI_high<0 and r.standardized_coefficient<0)
    candidates['priority_within_section_candidate']=candidates.pre_QC_priority_rule_passed&candidates.QC_masked_same_sign&candidates.QC_masked_interval_same_direction
    save(candidates,'candidate_relationships.tsv')

def main():
    p=argparse.ArgumentParser();p.add_argument('--step',choices=['audit','sensitivity','finalize','all'],default='all');args=p.parse_args()
    if args.step in ['audit','all']:audit()
    if args.step in ['sensitivity','all']:sensitivities()
    if args.step=='finalize':finalize_priorities()

if __name__=='__main__':main()
