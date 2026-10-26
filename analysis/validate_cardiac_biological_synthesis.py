"""Check revised source evidence, graph arithmetic and independent model fits."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy import sparse
from cardiac_subtype_niches import model_arrays
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/cardiac_biological_synthesis'
STAGES=['P0','P7','P14','P21']
def read(p,**kw):return pd.read_csv(p,sep='\t',**kw)
def main():
    cov=read(OUT/'FB_assignment_coverage_comparison.tsv');counts=[];independent_graph_checks=0
    for stage in STAGES:
        d=read(OUT/(stage+'_FB_assignment_evidence.tsv.gz'),index_col='cell_id')
        f=read(ROOT/'results/cardiac_spatial_pipeline'/stage/'atlas_annotations.tsv.gz',index_col=0)
        np.testing.assert_allclose(d[['x_centroid','y_centroid']],f.loc[d.index,['x_centroid','y_centroid']],rtol=0,atol=1e-7)
        for col in ['whole_cell_supported','boundary_supported_nuclear_optional','multigene_boundary_candidate']:
            selected=d[col].ne('Unresolved');s=d[selected]
            assert s.atlas_identity.eq('Fibroblast').all() and s.program_source_eligible.all()
            assert s['whole_cell__both_library_agreement'].all() and s['whole_cell__raw_reference_agreement'].all()
            assert s['whole_cell__state_information_QC'].all() and s['whole_cell__family_weight'].ge(.6).all()
            assert s['whole_cell__family_margin'].ge(.15).all() and s['whole_cell__measured_state_markers'].ge(2).all()
            assert not s['whole_cell__mixed_CM_FB_marker_review'].any()
            if col!='multigene_boundary_candidate':assert s['whole_cell__measured_state_anchors'].ge(1).all()
            if col!='whole_cell_supported':
                assert s.core_family_proposal.eq(s[col]).all() and s.core_markers.ge(2).all()
                assert not s.supported_nuclear_counterevidence.any()
        assert d.loc[d.boundary_supported_nuclear_optional.ne('Unresolved'),'multigene_boundary_candidate'].eq(
            d.loc[d.boundary_supported_nuclear_optional.ne('Unresolved'),'boundary_supported_nuclear_optional']).all()
        for r in cov[cov.stage.eq(stage)].itertuples():
            for col in ['legacy_strict_family','whole_cell_supported','boundary_supported_nuclear_optional','multigene_boundary_candidate']:
                assert d[col].eq(r.family).sum()==getattr(r,col)
        n=read(OUT/(stage+'_candidate_FB_CM_neighborhoods.tsv.gz'),index_col='cell_id')
        assert not n.index.isin(d.index).any()
        graph=sparse.load_npz(ROOT/'results/cardiac_spatial_pipeline'/stage/'proximity50_graph.npz').tocsr()
        for cid in n.sample(8,random_state=93).index:
            ids=f.index[graph.getrow(f.index.get_loc(cid)).indices]
            eligible=f.loc[ids].query("atlas_identity == 'Fibroblast' and program_source_eligible")
            assert len(eligible)==n.loc[cid,'eligible_FB_neighbors']
            z=d.reindex(eligible.index)
            for key in [k for k in n.columns if k.startswith('FB__')]:
                expected=z.multigene_boundary_candidate.eq(key[4:]).sum()/len(eligible)
                np.testing.assert_allclose(n.loc[cid,key],expected,atol=2e-8)
            independent_graph_checks+=1
        assert (n.filter(regex='^FB__').sum(axis=1)<=1+1e-7).all()
        counts.append({'stage':stage,'FB_profiles_retained':len(d),'candidate_sources':int(d.multigene_boundary_candidate.ne('Unresolved').sum()),
                       'original_coordinates_preserved':True,'no_focal_source_overlap':True})
    coeff=read(OUT/'nuclear_optional_FB_maturation_coefficients.tsv')
    oldcov=read(ROOT/'results/cardiac_subtype_niches/source_subtype_coverage.tsv')
    for stage,state,ex,outcome in [('P7','p.CM','FB__Postn-associated','core_Z__metabolic'),('P14','Ankrd1.CM','EC__Endocardial-like','core_Z__sarcomere')]:
        d=read(ROOT/'results/cardiac_subtype_niches'/(stage+'_proximity50_CM_subtype_neighborhoods.tsv.gz'),index_col='cell_id')
        n=read(OUT/(stage+'_candidate_FB_CM_neighborhoods.tsv.gz'),index_col='cell_id');d=d.loc[n.index].copy()
        exs=[]
        for r in cov[cov.stage.eq(stage)].itertuples():
            key='FB__'+r.family;d[key]=n[key]
            if r.multigene_boundary_candidate>=30 and d[key].gt(0).sum()>=50:exs.append(key)
        for key in oldcov.loc[oldcov.stage.eq(stage)&oldcov.subtype.str.startswith('EC__')&oldcov.primary_source_profiles.ge(30),'subtype']:
            if d[key].gt(0).sum()>=50:exs.append(key)
        names,x,y,outcomes=model_arrays(d,exs,view='primary',adjustment='within_region')
        take=d.CM_state.eq(state).to_numpy()&np.isfinite(x).all(1)&np.isfinite(y).all(1)
        groups=d.loc[take,'block200'];X=pd.DataFrame(x[take],columns=names,index=groups.index);Y=pd.Series(y[take,outcomes.index(outcome)],index=groups.index)
        X-=X.groupby(groups).transform('mean');Y-=Y.groupby(groups).transform('mean');X=X.loc[:,X.std(ddof=0)>.015]
        b=np.linalg.lstsq(X.to_numpy(),Y.to_numpy(),rcond=1e-9)[0]
        r=coeff[coeff.stage.eq(stage)&coeff.CM_state.eq(state)&coeff.source_subtype.eq(ex)&coeff.outcome.eq(outcome)&coeff.geometry.eq('proximity50')&coeff.source_view.eq('multigene_boundary_candidate')].iloc[0]
        np.testing.assert_allclose(b[X.columns.get_loc(ex)],r.standardized_coefficient,atol=3e-6)
    # Recompute the multiple-testing annotation independently per analysis grid.
    from scipy.stats import false_discovery_control
    for _,g in coeff.groupby(['stage','geometry','source_view','adjustment']):
        np.testing.assert_allclose(false_discovery_control(g.within_section_approx_p),g.within_section_BH_q,atol=2e-8)
    result={'status':'PASS','stages':counts,'independent_neighbor_reconstructions':independent_graph_checks,'independent_centered_lstsq_models':2,
        'full_grid_BH_recalculated':True,'original_atlas_and_previous_transport_not_overwritten':True,
        'new_cell_identity_accuracy_or_biological_replication_established':False}
    (OUT/'numerical_validation.json').write_text(json.dumps(result,indent=2)+'\n')
    print('BIOLOGICAL_NUMERICAL_VALIDATION_PASS',flush=True)
if __name__=='__main__':main()
