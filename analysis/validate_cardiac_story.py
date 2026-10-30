"""Independent numerical checks of new manuscript analyses."""
import json
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist,pdist
from scipy.optimize import linprog,linear_sum_assignment
from scipy import sparse
from scipy.stats import spearmanr
from cardiac_story_analysis import ROOT,OUT,STAGES,NICHES,RECON,UP,read,source,SOURCES

def main():
    obs=read(OUT/'OT_selected_observed_neighborhoods.tsv');pairs=read(OUT/'anatomy_stratified_neighborhood_pair_costs.tsv')
    plans=read(OUT/'adjacent_stage_neighborhood_correspondence_plans.tsv')
    assert len(obs)==96 and obs.cloud_index.nunique()==96 and len(pairs)==864
    meta=read(NICHES/'Wormhole_heldout_clouds.tsv');clouds={};checks=[]
    for stage in STAGES:
        print('STORY_VALIDATE',stage,flush=True)
        f=read(UP/stage/'atlas_annotations.tsv.gz',index_col=0);c=read(RECON/stage/'CM_outcomes_and_anatomy.tsv.gz',index_col=0)
        n=read(UP/stage/'niche_annotations.tsv.gz',index_col=0);sub=obs[obs.stage.eq(stage)]
        assert len(sub)==24
        z=f.loc[sub.cell_id];assert z.ventricular_CM_strict.all() and z.author_atlas_cell_state.str.contains('CM').all()
        assert c.loc[sub.cell_id,'reconstructed_Figure4e_CM'].all() and n.loc[sub.cell_id,'niche'].str.startswith('N2:').all()
        for stratum,g in sub.groupby('distance_stratum'):
            assert len(g)==12 and pdist(g[['x_um','y_um']]).min()>=100-1e-6
            near=f.loc[g.cell_id,'endocardial_proxy_distance_um'].le(100)
            assert (near if stratum.startswith('Endocardial') else ~near).all()
        with np.load(source(NICHES/stage/'nonCM_state_points.npz')) as b:x=b['points'];ids=b['focal_ids'].astype(str)
        member=sparse.load_npz(source(NICHES/stage/'state_cloud_membership.npz')).tocsr()
        for r in sub.itertuples():
            entry=meta.iloc[r.cloud_index];assert entry.stage==stage and entry.cell_id==r.cell_id
            assert ids[entry.focal_index]==r.cell_id
            ix=member.getrow(entry.focal_index).indices
            if len(ix)>32:ix=np.random.default_rng(947+int(entry.focal_index)).choice(ix,32,replace=False)
            clouds[r.cloud_index]=x[ix]
            assert len(ix)==r.sampled_source_profiles
        checks.append({'stage':stage,'observed_neighborhoods':len(sub),'eligibility_spacing_sources':'PASS'})
    np.testing.assert_allclose(pairs.exact_W2_squared,pairs.state_cost_on_optimal_plan+pairs.identity_cost_on_optimal_plan,rtol=2e-7)
    # Dense constraints independently rebuild one exact non-CM-cloud LP per
    # adjacent-stage×anatomy stratum (production uses sparse Kronecker matrices).
    independent_cloud_LPs=0;independent_neighborhood_assignments=0
    for keys,g in pairs.groupby(['stage_a','stage_b','distance_stratum']):
        r=g.iloc[17];x,y=clouds[int(r.cloud_a)],clouds[int(r.cloud_b)];na,nb=len(x),len(y)
        dense=np.zeros((na+nb,na*nb))
        for i in range(na):dense[i,i*nb:(i+1)*nb]=1
        for j in range(nb):dense[na+j,j::nb]=1
        lp=linprog(cdist(x,y,'sqeuclidean').ravel(),A_eq=dense,b_eq=np.r_[np.ones(na)/na,np.ones(nb)/nb],bounds=(0,None),method='highs')
        assert lp.success;np.testing.assert_allclose(lp.fun,r.exact_W2_squared,rtol=2e-7);independent_cloud_LPs+=1
        p=plans[plans.stage_a.eq(keys[0])&plans.stage_b.eq(keys[1])&plans.distance_stratum.eq(keys[2])]
        assert len(p)==144
        for name in ['balanced_correspondence_mass','entropic_correspondence_mass']:
            np.testing.assert_allclose(p.groupby('cloud_a')[name].sum(),1/12,atol=2e-8)
            np.testing.assert_allclose(p.groupby('cloud_b')[name].sum(),1/12,atol=2e-8)
        mat=g.pivot(index='cloud_a',columns='cloud_b',values='exact_W2_squared')
        ii,jj=linear_sum_assignment(mat.to_numpy())
        objective=mat.to_numpy()[ii,jj].mean()
        np.testing.assert_allclose((p.balanced_correspondence_mass*p.exact_cloud_distance).sum(),objective,rtol=2e-7)
        independent_neighborhood_assignments+=1
    metrics=read(OUT/'neighborhood_transport_metric_sensitivity.tsv')
    for r in metrics[metrics.metric.ne('primary')].itertuples():
        g=pairs[pairs.stage_a.eq(r.stage_a)&pairs.stage_b.eq(r.stage_b)&pairs.distance_stratum.eq(r.distance_stratum)]
        val=g[r.metric];np.testing.assert_allclose(spearmanr(g.exact_W2_squared,val).statistic,r.exact_cost_rank_correlation,atol=1e-8)
        exact=g.pivot(index='cloud_a',columns='cloud_b',values='exact_W2_squared')
        alt=g.pivot(index='cloud_a',columns='cloud_b',values=r.metric).reindex(index=exact.index,columns=exact.columns)
        agreement=np.mean(np.argmin(exact.to_numpy(),1)==np.argmin(alt.to_numpy(),1))
        np.testing.assert_allclose(agreement,r.nearest_exact_neighbor_agreement,atol=1e-8)
    with np.load(source(NICHES/'nonCM_state_PCA_model.npz')) as b:used=set(b['genes'].astype(str));excluded=set(b['excluded_genes'].astype(str))
    maturation=set(read('results/fibroblast_reassessment/CM_maturation_Figure4b_panel.tsv').gene)
    assert not used&maturation and maturation<=excluded
    fb=read(OUT/'fibroblast_library_matrix_programs.tsv');wide=fb[fb.level.eq('fibroblast_family')].pivot(index=['stage','library','group','program'],columns='view',values='mean_log1p_CPM')
    raw_corrected=spearmanr(wide.raw,wide.corrected).statistic
    assert raw_corrected>.95
    late=fb[fb.level.eq('fibroblast_family')&fb.program.eq('Basement matrix')&fb.stage.isin(['P14','P21'])&fb.group.isin(['Col8a1-associated','Gfpt2-associated'])]
    late=late.pivot(index=['stage','library','view'],columns='group',values='mean_log1p_CPM')
    assert (late['Col8a1-associated']>late['Gfpt2-associated']).all()
    candidates=read(OUT/'reviewed_fibroblast_candidate_state_continuities.tsv')
    assert candidates.loc[candidates.to_cluster.eq('P21_F0'),'to_family'].eq('Col8a1-associated').all()
    # Directly recalculate one full-RNA pseudobulk contrast using full gene
    # rows and actual matrix totals, independent of the selected-column writer.
    import anndata as ad
    ann=read('results/fibroblast_Xenium_reviewed_refit/RNA_clean_fibroblast_annotations.tsv.gz',index_col=0)
    a=ad.read_h5ad(source('/workspace/heart-atlas/results/cell_identities_subtypes/snrna_qc_annotated_subclustered.h5ad'),backed='r')
    panel='Col4a1 Col4a2 Col4a5 Col4a6 Lama4 Hspg2'.split();symbols=a.var.gene_symbol.astype(str).to_numpy()
    for family in ['Col8a1-associated','Gfpt2-associated']:
        selected=ann.index[ann.age.eq('P14')&ann['sample'].eq('GSM8827075_P14_1')&ann.fibroblast_family.eq(family)]
        ii=a.obs.index.get_indexer(selected);x=a.X[ii];pb=np.asarray(x.sum(0)).ravel();total=float(pb.sum())
        val=np.mean([np.log1p(pb[symbols==gene].sum()/total*1e6) for gene in panel])
        expected=fb[fb.level.eq('fibroblast_family')&fb.group.eq(family)&fb.stage.eq('P14')&fb.library.eq('GSM8827075_P14_1')&fb.view.eq('corrected')&fb.program.eq('Basement matrix')].mean_log1p_CPM.iloc[0]
        np.testing.assert_allclose(val,expected,rtol=2e-6)
    a.file.close()
    out={'status':'PASS','cohorts':checks,'exact_source_cloud_LPs_independently_recomputed':independent_cloud_LPs,
      'outer_neighborhood_LPs_checked_by_Hungarian_assignment':independent_neighborhood_assignments,
      'all_uniform_transport_marginals_checked':True,'entropic_marginals_checked':True,'all_cost_decompositions_checked':True,
      'rank_and_nearest_match_metrics_independently_recomputed':True,'CM_maturation_genes_withheld':True,
      'RNA_pseudobulks_independently_recomputed':2,'late_basement_contrast_supported_in_each_library_and_raw_corrected_view':True,
      'RNA_raw_corrected_program_Spearman':raw_corrected,'P21_F0_current_Col8a1_family_preserved':True,
      'biological_fates_validated':False,'new_biological_replicates':0,'scope':'Numerical/data provenance; no anatomical, causal or lineage validation'}
    (OUT/'numerical_validation.json').write_text(json.dumps(out,indent=2)+'\n')
    save_checks=pd.DataFrame(checks);save_checks.to_csv(OUT/'cohort_validation.tsv',sep='\t',index=False)
    p=OUT/'analysis_source_files.json';SOURCES.update(json.loads(p.read_text()));p.write_text(json.dumps(sorted(SOURCES),indent=2)+'\n')
    print('STORY_NUMERICAL_VALIDATION_PASS',flush=True)

if __name__=='__main__':main()
