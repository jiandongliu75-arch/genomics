"""Independent distance, regression, multiplicity and matched-score checks."""
from pathlib import Path
import json,hashlib
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import false_discovery_control,t
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/capillary_maturation_retest'
def main():
 p=pd.read_csv(OUT/'proximity_regressions.tsv',sep='\t')
 m=pd.read_csv(OUT/'PTPRM_and_Vegfa_regressions.tsv',sep='\t')
 assert p.primary_analysis.sum()==8 and m.primary_PTPRM_analysis.sum()==60
 for d,mask,q in [(p,p.primary_analysis,'primary_BH_q'),(m,m.primary_PTPRM_analysis,'PTPRM_BH_q')]:
  np.testing.assert_allclose(false_discovery_control(d.loc[mask,'within_section_p']),d.loc[mask,q],rtol=1e-6,atol=1e-12)
 distances=[];coefs=[]
 for stage in ['P0','P7','P14','P21']:
  c=pd.read_csv(OUT/(stage+'_CM_retest_values.tsv.gz'),sep='\t',index_col=0)
  assert c.index.is_unique and np.isfinite(c.core_maturation21).all()
  f=pd.read_csv(ROOT/f'results/cardiac_spatial_pipeline/{stage}/atlas_annotations.tsv.gz',sep='\t',index_col=0)
  mask=f.atlas_identity.eq('Endothelial')&f.program_source_eligible&f.EC_nuclear_type.eq('Capillary-like')
  distance=cKDTree(f.loc[mask,['x_centroid','y_centroid']]).query(c[['x_centroid','y_centroid']])[0]
  error=float(np.max(abs(distance-c.distance_QC_capillary)));assert error<2e-4
  distances.append(dict(stage=stage,maximum_distance_error_um=error,CMs=len(c)))
  sub=c[c.strict];n=len(sub)
  # Reconstruct the full design without production QR or fitting helpers.
  numeric=np.column_stack([np.log1p(sub.capture),np.log1p(sub.cell_area),sub.nucleus_missing.astype(float),
   np.log1p(sub.distance_to_large_void_um),np.log1p(sub.distance_to_outer_um)])
  dummy=pd.get_dummies(sub[['CM_state','morphology_region']],dtype=float).to_numpy()
  exposure=sub.near_QC_capillary.to_numpy();exposure=(exposure-exposure.mean())/exposure.std()
  y=sub.core_maturation21.to_numpy();y=(y-y.mean())/y.std()
  X=np.column_stack([np.ones(n),numeric,dummy,exposure]);beta=np.linalg.lstsq(X,y,rcond=1e-10)[0]
  target=p[p.primary_analysis&p.stage.eq(stage)&p.exposure.eq('near_QC_capillary')].iloc[0]
  error=abs(beta[-1]-target.beta_SD);assert error<1e-7
  # Independent partial-regression clustered SE using the FWL residual.
  C=X[:,:-1];xr=exposure-C@np.linalg.lstsq(C,exposure,rcond=1e-10)[0]
  residual=y-X@beta;block=(sub.x_centroid//400).astype(int).astype(str)+'_'+(sub.y_centroid//400).astype(int).astype(str)
  score=pd.Series(xr*residual,index=sub.index).groupby(block).sum().to_numpy();g=len(score)
  se=np.sqrt(np.sum(score**2)/(xr@xr)**2*g/(g-1)*(n-1)/(n-int(target.design_rank)))
  ci=target.beta_SD+t.ppf(.975,g-1)*np.array([-se,se])
  np.testing.assert_allclose(ci,[target.CI_low,target.CI_high],atol=1e-7)
  coefs.append(dict(stage=stage,independent_slope_error=float(error),spatial_blocks=g,FWL_cluster_interval_verified=True))
 old=pd.read_csv(ROOT/'results/CM_maturation21_nonCM_context/matched_capillary_maturation_results.tsv',sep='\t')
 new=pd.read_csv(OUT/'independently_recomputed_matched_differences.tsv',sep='\t')
 errors=[]
 for r in new.itertuples():
  a=old[old.stage.eq(r.stage)&old.cohort.eq(r.cohort)&old.outcome.eq(r.view+'_Z__maturation21')]
  assert len(a)==1
  errors.append(abs(a.iloc[0].raw_score_difference-r.raw_score_difference))
 assert max(errors)<1e-7
 checks=pd.read_csv(OUT/'PTPRM_count_reconstruction_checks.tsv',sep='\t');assert len(checks)==8 and checks.status.eq('PASS').all()
 s=pd.read_csv(OUT/'flexible_capture_sensitivity.tsv',sep='\t')
 a=s[s.stage.eq('P7')&s.exposure.eq('near_QC_capillary')&s.adjustment.eq('within_region_state_spline')]
 assert a.CI_low.gt(0).all()
 loo=p[p.stage.eq('P7')&p.outcome.str.contains('_omit_')];assert len(loo)==84 and loo.CI_low.gt(0).all()
 docs=json.loads((OUT/'document_validation.json').read_text());assert docs['status']=='PASS'
 methods=json.loads((OUT/'methods.json').read_text());assert methods['primary_proximity_tests']==8 and methods['primary_PTPRM_tests']==60
 result=dict(status='PASS',primary_proximity_tests=8,primary_PTPRM_tests=60,independent_distances=distances,
  independent_OLS_and_FWL_intervals=coefs,matched_contrasts_checked=len(errors),maximum_matched_difference_error=max(errors),
  count_rebuilt_PTPRM_models=8,P7_positive_leave_one_out_relationships=84,P7_flexible_capture_and_depth_tests_positive=True,
  no_new_biological_replicates=True,documents=docs)
 (OUT/'validation.json').write_text(json.dumps(result,indent=2)+'\n')
 files=sorted(p for p in OUT.iterdir() if p.is_file() and p.name!='file_manifest.tsv')
 pd.DataFrame([dict(path=p.name,bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in files]).to_csv(OUT/'file_manifest.tsv',sep='\t',index=False)
 print(json.dumps(result,indent=2))
if __name__=='__main__':main()
