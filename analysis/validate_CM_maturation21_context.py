"""Independent neighbor counts, regression reconstruction and matching checks."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import sparse
import assess_CM_maturation21_neighborhoods as m

def main():
    independent=[]; neighbor_checks=0;ols_checks=0
    coefficients=m.read(m.OUT/'nonCM_maturation_associations.tsv')
    for stage in m.STAGES:
        f,d,coverage=m.sources_and_base(stage,'proximity50')
        saved=m.read(m.OUT/stage/'CM_nonCM_neighborhoods.tsv.gz',index_col='cell_id')
        fb=m.read(f'results/fibroblast_multimethod_identity/{stage}_multimethod_FB_annotations.tsv.gz',index_col='cell_id')
        g=sparse.load_npz(m.ROOT/f'results/cardiac_spatial_pipeline/{stage}/proximity50_graph.npz').tocsr()
        assert np.count_nonzero(g.diagonal())==0
        rng=np.random.default_rng(m.SEED+int(stage[1:]));sample=rng.choice(d.index,size=32,replace=False)
        for cid in sample:
            i=f.index.get_loc(cid);neighbors=f.iloc[g.indices[g.indptr[i]:g.indptr[i+1]]]
            eligible=neighbors[neighbors.program_source_eligible]
            resolved=int(neighbors.atlas_resolved.sum())
            for source in coverage.source:
                kind,name=source.split('__');den=resolved
                if kind=='Broad':num=int(eligible.atlas_identity.eq(name).sum())
                elif kind=='EC':
                    e=eligible[eligible.atlas_identity.eq('Endothelial')];den=len(e);num=int(e.EC_nuclear_type.eq(name).sum())
                else:
                    e=eligible[eligible.atlas_identity.eq('Fibroblast')];den=len(e)
                    ids=e.index.intersection(fb.index);num=int(fb.loc[ids,'recommended_candidate'].eq(name+'-associated').sum())
                expected=num/max(den,1);actual=saved.loc[cid,source]
                assert abs(expected-actual)<1e-8,(stage,cid,source,expected,actual)
                neighbor_checks+=1
        z=d[d.analysis_eligible].copy()
        candidates=coefficients[coefficients.stage.eq(stage)&coefficients.primary_analysis]
        chosen=candidates[candidates.source.isin(['Broad__Fibroblast','Broad__Smooth muscle','EC__Capillary-like'])]
        chosen=chosen[chosen.outcome.eq('core_Z__maturation21') | (chosen.source.eq('Broad__Fibroblast')&chosen.outcome.eq('core_Z__metabolic'))]
        for _,row in chosen.iterrows():
            x=m.nuisance(z,row.source).astype(float);e=z[row.source].astype(float);y=z[row.outcome].astype(float)
            finite=np.isfinite(x).all(axis=1)&np.isfinite(e)&np.isfinite(y)
            x=(x-x.mean())/np.maximum(x.std(ddof=0),1e-6)
            e=(e-e.mean())/max(e.std(ddof=0),1e-6);y=(y-y.mean())/max(y.std(ddof=0),1e-6)
            x=x.loc[finite];e=e.loc[finite];y=y.loc[finite]
            group=z.loc[finite,'block200']+'::'+z.loc[finite,'CM_state']
            x=x-x.groupby(group).transform('mean');e=e-e.groupby(group).transform('mean');y=y-y.groupby(group).transform('mean')
            keep=x.std(ddof=0)>.015
            design=np.column_stack([x.loc[:,keep],e]);b=np.linalg.lstsq(design,y,rcond=1e-8)[0][-1]
            error=abs(b-row.coefficient);assert error<2e-7,(stage,row.source,error)
            independent.append({'stage':stage,'source':row.source,'outcome':row.outcome,'maximum_absolute_error':error,'status':'PASS'})
            ols_checks+=1
    pairs=m.read(m.OUT/'QC_capillary_refitted_pairs.tsv.gz')
    for stage,p in pairs.groupby('stage'):
        assert not p.high_apposition_CM.duplicated().any() and not p.zero_apposition_CM.duplicated().any()
        assert set(p.high_apposition_CM).isdisjoint(p.zero_apposition_CM)
        assert p.centroid_separation_um.max()<=100+1e-6 and p.matching_distance.max()<=4+1e-6
        assert not p.matching_uses_maturation_score.any()
        c=m.read(f'results/cardiac_spatial_reconstruction/{stage}/CM_outcomes_and_anatomy.tsv.gz',index_col=0)
        for field in ['author_atlas_cell_state','anatomy_block_200um']:
            assert np.array_equal(c.loc[p.high_apposition_CM,field].to_numpy(),c.loc[p.zero_apposition_CM,field].to_numpy())
    m.save(pd.DataFrame(independent),'independent_regression_checks.tsv')
    result={'status':'PASS','independent_neighbor_fraction_checks':neighbor_checks,
            'independently_reconstructed_regression_coefficients':ols_checks,'QC_capillary_pairs_checked':len(pairs),
            'QC_matching_uses_scores':False,'QC_pairs_same_CM_state_and_local_region':True}
    (m.OUT/'independent_validation.json').write_text(json.dumps(result,indent=2)+'\n')
    print('INDEPENDENT_MATURATION21_VALIDATION_PASS',json.dumps(result),flush=True)

if __name__=='__main__':main()
