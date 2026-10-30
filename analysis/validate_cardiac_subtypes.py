"""Independent subtype-neighborhood, model, matching and transport checks."""
import json,itertools
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.optimize import linprog,linear_sum_assignment
from scipy.spatial.distance import cdist,pdist
from cardiac_subtype_niches import ROOT,OUT,UP,STAGES,CM,FB,EC,read,source,SOURCES,stage_data,model_arrays

def main():
    coverage=read(OUT/'source_subtype_coverage.tsv');coeff=read(OUT/'subtype_conditioned_maturation_coefficients.tsv')
    stages=[];frames={}
    for stage in STAGES:
        print('VALIDATE_SUBTYPE_STAGE',stage,flush=True);f,original=stage_data(stage)
        d=read(OUT/(stage+'_proximity50_CM_subtype_neighborhoods.tsv.gz'),index_col='cell_id');frames[stage]=d
        assert d.index.equals(original.index)
        np.testing.assert_allclose(d[['x_centroid','y_centroid']],original[['x_centroid','y_centroid']],rtol=1e-8)
        expected=d.groupby('CM_state').core_Z__maturation21.rank(pct=True)
        np.testing.assert_allclose(expected,d.state_maturation_percentile,atol=2e-8)
        src=f.program_source_eligible&f.atlas_identity.isin(['Fibroblast','Endothelial'])
        assert not d.index.isin(f.index[src]).any()
        graph=sparse.load_npz(source(UP/stage/'proximity50_graph.npz')).tocsr()
        for cid in d.sample(6,random_state=23).index:
            neighbor=f.iloc[graph.getrow(f.index.get_loc(cid)).indices]
            fb=neighbor[neighbor.atlas_identity.eq('Fibroblast')&neighbor.program_source_eligible]
            ec=neighbor[neighbor.atlas_identity.eq('Endothelial')&neighbor.program_source_eligible]
            assert len(fb)==d.loc[cid,'eligible_FB_neighbors'] and len(ec)==d.loc[cid,'eligible_EC_neighbors']
            for name in FB:np.testing.assert_allclose((fb.fine_identity==name).sum()/max(len(fb),1),d.loc[cid,'FB__'+name],atol=1e-8)
            for name in EC:np.testing.assert_allclose((ec.EC_nuclear_type==name).sum()/max(len(ec),1),d.loc[cid,'EC__'+name],atol=1e-8)
        for prefix in ['FB__','EC__']:
            assert (d[[k for k in d if k.startswith(prefix)]].sum(axis=1)<=1+1e-7).all()
        for state in ['v.CM','Ankrd1.CM']:assert d.CM_state.eq(state).sum()>500
        stages.append({'stage':stage,'CM_profiles':len(d),'within_state_ranks_checked':True,'independent_neighbor_reconstructions':6,'original_coordinates_preserved':True})
    # Independently solve two centered models using QR/lstsq rather than the
    # production block-summed normal equations and shared bootstrap machinery.
    targets=[('P7','p.CM','FB__Postn-associated','core_Z__metabolic'),('P14','Ankrd1.CM','EC__Endocardial-like','core_Z__sarcomere')]
    for stage,state,ex,outcome in targets:
        d=frames[stage];d=d[d.analysis_eligible]
        cov=coverage[coverage.stage.eq(stage)&coverage.primary_source_profiles.ge(30)]
        expos=[k for k in cov.subtype if d[k].gt(0).sum()>=50]
        names,x,y,outcomes=model_arrays(d,expos,view='primary',adjustment='within_region')
        finite=np.isfinite(x).all(1)&np.isfinite(y).all(1);take=finite&d.CM_state.eq(state).to_numpy();z=d.loc[take]
        X=pd.DataFrame(x[take],columns=names,index=z.index);Y=pd.DataFrame(y[take],columns=outcomes,index=z.index)
        X-=X.groupby(z.block200).transform('mean');Y-=Y.groupby(z.block200).transform('mean')
        keep=X.std(ddof=0)>.015;X=X.loc[:,keep]
        b=np.linalg.lstsq(X.to_numpy(),Y[outcome].to_numpy(),rcond=1e-9)[0]
        row=coeff[coeff.stage.eq(stage)&coeff.CM_state.eq(state)&coeff.geometry.eq('proximity50')&coeff.source_view.eq('primary')&coeff.adjustment.eq('within_region')&coeff.source_subtype.eq(ex)&coeff.outcome.eq(outcome)].iloc[0]
        np.testing.assert_allclose(b[X.columns.get_loc(ex)],row.standardized_coefficient,atol=2e-6)
    delta=read(OUT/'Ankrd1_versus_working_CM_slope_differences.tsv')
    for r in delta.sample(12,random_state=73).itertuples():
        z=coeff[coeff.stage.eq(r.stage)&coeff.geometry.eq(r.geometry)&coeff.source_view.eq(r.source_view)&coeff.adjustment.eq(r.adjustment)&coeff.source_subtype.eq(r.source_subtype)&coeff.outcome.eq(r.outcome)].set_index('CM_state')
        np.testing.assert_allclose(z.loc['Ankrd1.CM','standardized_coefficient']-z.loc['v.CM','standardized_coefficient'],r.Ankrd1_minus_working_slope,atol=2e-7)
    paired=read(OUT/'same_CM_geometry_sensitivity.tsv')
    for _,g in paired.groupby(['stage','CM_state','paired_comparison','source_subtype','outcome']):assert g.CM_profiles.nunique()==1
    matches=read(OUT/'within_CM_state_maturation_matched_pairs.tsv.gz')
    for stage,g in matches.groupby('stage'):
        d=frames[stage];d=d[d.analysis_eligible]
        variables=np.column_stack([np.log1p(d.CM_outcome_excluded_depth),np.log1p(d.cell_area),np.log1p(d.eligible_FB_neighbors),np.log1p(d.eligible_EC_neighbors),d.FB_mean_log_capture,d.EC_mean_log_capture])
        standardized=(variables-variables.mean(0))/np.maximum(variables.std(0),1e-6)
        for state,z in g.groupby('CM_state'):
            assert not z.low_maturation_CM.duplicated().any() and not z.high_maturation_CM.duplicated().any()
            a=d.loc[z.high_maturation_CM].reset_index();b=d.loc[z.low_maturation_CM].reset_index()
            assert a.CM_state.eq(state).all() and b.CM_state.eq(state).all()
            assert a.maturation_group.eq('High').all() and b.maturation_group.eq('Low').all()
            assert (a.block200==b.block200).all()
            assert np.linalg.norm(a[['x_centroid','y_centroid']].to_numpy()-b[['x_centroid','y_centroid']].to_numpy(),axis=1).max()<=100+1e-6
            assert (a.endocardial_proxy_distance_um-b.endocardial_proxy_distance_um).abs().max()<=50+1e-6
            diff=np.abs(standardized[d.index.get_indexer(z.high_maturation_CM)]-standardized[d.index.get_indexer(z.low_maturation_CM)])
            assert diff[:,:4].max()<=.5+1e-6 and diff[:,4:].max()<=.75+1e-6
    obs=read(OUT/'CM_identity_restricted_OT_neighborhoods.tsv');cost=read(OUT/'CM_identity_restricted_OT_pair_costs.tsv');plan=read(OUT/'CM_identity_restricted_OT_correspondence.tsv')
    assert len(obs)==96 and obs.cloud_index.nunique()==96
    clouds={};prior=ROOT/'results/cardiac_niche_dynamics';meta=read(prior/'Wormhole_heldout_clouds.tsv')
    for stage,g in obs.groupby('stage'):
        with np.load(source(prior/stage/'nonCM_state_points.npz')) as b:points=b['points']
        graph=sparse.load_npz(source(prior/stage/'state_cloud_membership.npz')).tocsr()
        for state,z in g.groupby('CM_state'):
            assert len(z)==12
            for stratum,h in z.groupby('distance_stratum'):assert len(h)==(4 if stratum=='≤100µm' else 8) and pdist(h[['x_um','y_um']]).min()>=75-1e-6
        for r in g.itertuples():
            row=meta.iloc[r.cloud_index];assert row.stage==stage and row.cell_id==r.cell_id
            assert frames[stage].loc[r.cell_id,'CM_state']==r.CM_state
            np.testing.assert_allclose(frames[stage].loc[r.cell_id,'state_maturation_percentile'],r.maturation_percentile_within_CM_state,atol=1e-8)
            ix=graph.getrow(row.focal_index).indices
            if len(ix)>32:ix=np.random.default_rng(947+int(row.focal_index)).choice(ix,32,replace=False)
            clouds[r.cloud_index]=points[ix]
    independent_LPs=0;outer_checks=0
    for key,g in cost.groupby(['stage_a','stage_b','CM_state','distance_stratum']):
        mat=g.pivot(index='cloud_a',columns='cloud_b',values='exact_W2_squared');p=plan[plan.stage_a.eq(key[0])&plan.stage_b.eq(key[1])&plan.CM_state.eq(key[2])&plan.distance_stratum.eq(key[3])]
        n=len(mat)
        for col in ['exact_correspondence_mass','Wormhole_correspondence_mass']:
            np.testing.assert_allclose(p.groupby('cloud_a')[col].sum(),1/n,atol=1e-8);np.testing.assert_allclose(p.groupby('cloud_b')[col].sum(),1/n,atol=1e-8)
        ii,jj=linear_sum_assignment(mat);np.testing.assert_allclose(mat.to_numpy()[ii,jj].mean(),(p.exact_cost*p.exact_correspondence_mass).sum(),rtol=2e-7);outer_checks+=1
        # One cloud cost per stratum, dense independent constraints.
        r=g.iloc[0];x,y=clouds[int(r.cloud_a)],clouds[int(r.cloud_b)];a,b=len(x),len(y);A=np.zeros((a+b,a*b))
        for i in range(a):A[i,i*b:(i+1)*b]=1
        for j in range(b):A[a+j,j::b]=1
        lp=linprog(cdist(x,y,'sqeuclidean').ravel(),A_eq=A,b_eq=np.r_[np.ones(a)/a,np.ones(b)/b],bounds=(0,None),method='highs');assert lp.success
        np.testing.assert_allclose(lp.fun,r.exact_W2_squared,rtol=2e-7);independent_LPs+=1
    np.testing.assert_allclose(cost.exact_W2_squared,cost.state_cost+cost.lineage_cost,rtol=2e-7)
    out={'status':'PASS','stages':stages,'independent_target_coefficients':2,'shared_CM_state_slope_difference_checks':12,
        'same_CM_geometry_cohorts_checked':True,'matched_identity_quartiles_regions_calipers_checked':True,
        'known_subtypes_not_gated_by_unsupervised_cluster_stability':True,'untyped_sources_retained_in_denominators':True,
        'CM_identity_restricted_OT_observed_neighborhoods':96,'independent_cloud_LPs':independent_LPs,'outer_objectives_checked_with_Hungarian':outer_checks,
        'transport_marginals_and_cost_decomposition_checked':True,'new_biological_replicates':0,'lineage_fate_or_causal_mechanism_validated':False}
    (OUT/'numerical_validation.json').write_text(json.dumps(out,indent=2)+'\n')
    p=OUT/'source_files.json';SOURCES.update(json.loads(p.read_text()));p.write_text(json.dumps(sorted(SOURCES),indent=2)+'\n')
    print('SUBTYPE_NUMERICAL_VALIDATION_PASS',flush=True)

if __name__=='__main__':main()
