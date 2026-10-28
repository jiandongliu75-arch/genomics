"""Check geometry, match eligibility and independently refit three primary slopes."""
from pathlib import Path
import json,hashlib
import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control
from retest_capillary_maturation import fit
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/capillary_association_robustness'
def main():
 r=pd.read_csv(OUT/'P7_robustness_regressions.tsv',sep='\t');primary=r[r.primary_analysis]
 assert len(primary)==3
 np.testing.assert_allclose(false_discovery_control(primary.wild_cluster_null_p),primary.primary_BH_q,atol=1e-12)
 d=pd.read_csv(ROOT/'results/capillary_maturation_retest/P7_CM_retest_values.tsv.gz',sep='\t',index_col=0);d=d[d.strict]
 g=pd.read_csv(OUT/'P7_added_geometry_values.tsv.gz',sep='\t',index_col=0);assert d.index.equals(g.index)
 d['near_boundary']=-np.log1p(g.boundary_um);d['near_tissue_path']=-np.log1p(g.path_um)
 reachable=np.isfinite(g.path_um)
 assert (g.loc[reachable,'path_um']+2e-4>=d.loc[reachable,'distance_QC_capillary']).all()
 checks=[]
 for row in primary.itertuples():
  independent=fit(d.loc[np.isfinite(d[row.exposure])],[row.exposure],['core_maturation21'],'within_region_state_spline')[0]
  err=abs(independent['beta_SD']-row.beta_SD);assert err<1e-7
  np.testing.assert_allclose([independent['CI_low'],independent['CI_high']],[row.CI_low,row.CI_high],atol=1e-7)
  checks.append(dict(exposure=row.exposure,independent_coefficient_error=err,independent_interval_verified=True))
 grid=r[r.analysis.eq('grid_sensitivity')];assert len(grid)==36 and grid.beta_SD.gt(0).all() and grid.CI_low.gt(0).all()
 pp=pd.read_csv(OUT/'P7_local_near_far_pairs.tsv.gz',sep='\t')
 assert len(pp)==916 and pp.near_CM.is_unique and pp.far_CM.is_unique and set(pp.near_CM).isdisjoint(pp.far_CM)
 a=d.loc[pp.near_CM];b=d.loc[pp.far_CM]
 assert a.CM_state.to_list()==b.CM_state.to_list() and a.morphology_region.to_list()==b.morphology_region.to_list()
 assert a.distance_QC_capillary.le(25).all() and b.distance_QC_capillary.ge(50).all()
 physical=np.linalg.norm(a[['x_centroid','y_centroid']].to_numpy()-b[['x_centroid','y_centroid']].to_numpy(),axis=1)
 assert physical.max()<=100+1e-6
 np.testing.assert_allclose(physical,pp.physical_pair_distance_um,atol=1e-5)
 xyA=np.floor(a[['x_centroid','y_centroid']].to_numpy()/200);xyB=np.floor(b[['x_centroid','y_centroid']].to_numpy()/200);assert np.array_equal(xyA,xyB)
 balance=[]
 for col in ['capture','cell_area']:
  full=np.log1p(d[col]);delta=(np.log1p(a[col]).to_numpy()-np.log1p(b[col]).to_numpy())/full.std(ddof=0)
  assert abs(delta).max()<=.25+1e-6
  balance.append(dict(variable=col,mean_standardized_pair_difference=delta.mean(),maximum_absolute_pair_difference=abs(delta).max()))
 matched=pd.read_csv(OUT/'P7_local_matched_outcomes.tsv',sep='\t')
 for row in matched.itertuples():
  values=a[row.outcome].to_numpy()-b[row.outcome].to_numpy();values=values[np.isfinite(values)]
  point=values.mean()/d[row.outcome].std(ddof=0)
  assert len(values)==row.matched_pairs and abs(point-row.standardized_difference)<1e-8
 assert matched.CI_low.lt(0).all() and matched.CI_high.gt(0).all()
 pd.DataFrame(balance).to_csv(OUT/'P7_matching_balance.tsv',sep='\t',index=False)
 docs=json.loads((OUT/'document_validation.json').read_text());assert docs['status']=='PASS'
 inputs=pd.read_csv(OUT/'input_checksums.tsv',sep='\t')
 for row in inputs.itertuples():assert hashlib.sha256((ROOT/row.path).read_bytes()).hexdigest()==row.sha256
 result=dict(status='PASS',primary_refits=checks,geometry_triangle_check=True,grid_tests=36,eligible_disjoint_matched_pairs=916,
  all_matched_outcome_points_independently_recomputed=True,matching_balance=balance,documents=docs,
  input_files_verified=len(inputs),new_biological_replicates=0)
 (OUT/'validation.json').write_text(json.dumps(result,indent=2)+'\n')
 files=sorted(p for p in OUT.iterdir() if p.is_file() and p.name!='file_manifest.tsv')
 pd.DataFrame([dict(path=p.name,bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in files]).to_csv(OUT/'file_manifest.tsv',sep='\t',index=False)
 print(json.dumps(result,indent=2))
if __name__=='__main__':main()
