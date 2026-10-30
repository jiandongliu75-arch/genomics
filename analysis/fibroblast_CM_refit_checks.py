"""Matched-radius tests, spatial descriptions and independent numerical checks."""
import json
import numpy as np
import pandas as pd
from scipy import sparse
from refit_fibroblast_CM_relationships import *

def loaded(stage,geom):
 d=model_base(stage,geom)
 with np.load(OUT/(stage+'_'+geom+'_neighborhood_exposures.npz')) as b:
  assert np.array_equal(b['cell_ids'],d.index.to_numpy(str));e=pd.DataFrame({k:b[k] for k in b.files if k!='cell_ids'},index=d.index)
 return pd.concat([d.drop(columns=d.columns.intersection(e.columns)),e],axis=1)
def eligible(d):return d.analysis_eligible&(d.refined__informative_neighbors>=1)&(d.core__informative_neighbors>=1)&(d.native__informative_neighbors>=1)
def paired():
 rows=[];diag=[];coverage=read(OUT/'fibroblast_assignment_coverage.tsv')
 for stage in STAGES:
  a=loaded(stage,'proximity50')
  for other in ['proximity30','contact']:
   b=loaded(stage,other);ids=a.index[eligible(a)&eligible(b)]
   targets=[('identity','candidate','candidate__'+f) for f in coverage.loc[coverage.stage.eq(stage)&coverage.label_layer.eq('candidate')&coverage.candidate_profiles.ge(30),'family']]
   targets += [('program','refined','refined__'+m) for m in MODULES]
   for geom,d in [('proximity50',a.loc[ids]),(other,b.loc[ids])]:
    for kind,layer,target in targets:
     r,_,di=fit_target(d,target,kind,layer,'within_region',geom)
     rows.extend([{'paired_geometry':other,'same_CM_profiles':True,**x} for x in r]);diag.extend([{'paired_geometry':other,**x} for x in di])
   print('REFIT_PAIRED',stage,other,len(ids),flush=True)
 save(pd.DataFrame(rows),'same_CM_radius_and_boundary_sensitivity.tsv');save(pd.DataFrame(diag),'same_CM_sensitivity_designs.tsv');export_inputs()

def spatial():
 from name_fibroblasts_and_describe_neighbors import matched_context
 rows=[];pairs=[];balances=[];profiles=read(OUT/'fibroblast_per_profile_programs.tsv.gz');summary=[]
 for stage in STAGES:
  c=read(f'results/fibroblast_names_and_neighborhoods/{stage}/FB_physical_neighbors.tsv.gz',index_col=0)
  labels=read(f'results/fibroblast_multimethod_identity/{stage}_multimethod_FB_annotations.tsv.gz',index_col='cell_id')
  for layer,col in [('candidate','recommended_candidate'),('crosschecked','crosschecked_candidate')]:
   cc=c.copy();cc['high_confidence_family']=labels[col].where(labels.program_source_eligible,'Unresolved')
   if stage in ['P7','P14']:
    r,p,b=matched_context(stage,cc)
    rows.extend([{'new_label_layer':layer,**x} for x in r]);pairs.extend([{'new_label_layer':layer,**x} for x in p]);balances.extend([{'new_label_layer':layer,**x} for x in b])
   for fam,g in cc[cc.high_confidence_family.ne('Unresolved')].groupby('high_confidence_family'):
    for metric in ['arterial_within_50um','arterial_distance_um','fraction_Cardiomyocyte_50um','fraction_Endothelial_50um','fraction_Pericyte_50um','fraction_Smooth muscle_50um','fraction_Immune_50um']:
     if metric in g:summary.append({'stage':stage,'label_layer':layer,'family':fam,'measure':metric,'profiles':len(g),'mean':g[metric].mean(),'median':g[metric].median()})
  p=profiles[profiles.stage.eq(stage)&profiles.view.eq('refined')].set_index('cell_id')
  for fam,g in labels[labels.program_source_eligible].groupby('recommended_candidate'):
   for mod in MODULES:summary.append({'stage':stage,'label_layer':'candidate','family':fam,'measure':'program__'+mod,'profiles':len(g),'mean':p.loc[g.index,'program__'+mod].mean(),'median':p.loc[g.index,'program__'+mod].median()})
 save(pd.DataFrame(rows),'updated_matched_fibroblast_spatial_context.tsv');save(pd.DataFrame(pairs),'updated_matched_fibroblast_pairs.tsv.gz');save(pd.DataFrame(balances),'updated_spatial_matching_balance.tsv');save(pd.DataFrame(summary),'updated_fibroblast_spatial_and_program_summary.tsv');export_inputs()

def validation():
 labels=read(OUT/'fibroblast_assignment_coverage.tsv');programs=read(OUT/'fibroblast_per_profile_programs.tsv.gz')
 coefficients=read(OUT/'CM_fibroblast_refitted_coefficients.tsv');diagnostics=read(OUT/'model_design_diagnostics.tsv')
 counts=read('results/fibroblast_multimethod_identity/final_candidate_coverage.tsv')
 for _,r in labels[labels.label_layer.eq('candidate')&labels.reference_assessable].iterrows():
  oldrow=counts[counts.stage.eq(r.stage)&counts.family.eq(r.family+'-associated')].iloc[0];assert int(r.candidate_profiles)==int(oldrow.candidate_profiles)
 assert labels.loc[~labels.reference_assessable,'candidate_profiles'].isna().all()
 assert not (set(sum(PANELS.values(),[]))&set(read('results/fibroblast_reassessment/CM_maturation_Figure4b_panel.tsv').gene))
 # Verify exact graph/ID aggregation independently, without sparse multiplication.
 checked=0;rng=np.random.default_rng(947)
 for stage in STAGES:
  f=read(f'results/cardiac_spatial_pipeline/{stage}/atlas_annotations.tsv.gz',index_col=0)
  lab=read(f'results/fibroblast_multimethod_identity/{stage}_multimethod_FB_annotations.tsv.gz',index_col='cell_id')
  for geom in GEOMETRIES:
   path=ROOT/f'results/cardiac_spatial_pipeline/{stage}/{geom}_graph.npz';g=sparse.load_npz(path).tocsr()
   with np.load(OUT/(stage+'_'+geom+'_neighborhood_exposures.npz')) as b:
    for i in rng.choice(len(b['cell_ids']),size=24,replace=False):
     cid=str(b['cell_ids'][i]);ii=f.index.get_loc(cid);neighbors=f.index[g.getrow(ii).indices];fb=neighbors.intersection(lab.index)
     good=fb[lab.loc[fb,'program_source_eligible'].to_numpy(bool)]
     assert len(good)==b['eligible_FB_neighbors'][i]
     for fam in FAMILIES:
      val=lab.loc[good,'recommended_candidate'].eq(fam+'-associated').sum()/max(len(good),1)
      assert abs(val-b['candidate__'+fam][i])<1e-9
     for view in ['refined','core','native','depthmatched']:
      p=programs[programs.stage.eq(stage)&programs.view.eq(view)].set_index('cell_id');use=good[p.loc[good,'informative'].to_numpy(bool)]
      for mod in MODULES:
       val=p.loc[use,'program__'+mod].mean() if len(use) else 0
       assert abs(val-b[view+'__'+mod][i])<1e-8
       checked+=1
 # Independently residualize outcome/target against nuisance and block intercepts.
 regressions=0
 for stage,kind,layer,target,state in [('P7','identity','candidate','candidate__Postn','p.CM'),('P14','program','refined','refined__Basement_matrix','Ankrd1.CM'),('P21','program','refined','refined__Gfpt2.Fib','v.CM')]:
  d=loaded(stage,'proximity50');d=d[eligible(d)].copy();cov=nuisance(d,'within_region',layer if kind=='program' else None)
  x=cov.to_numpy(float);x=(x-np.nanmean(x,0))/np.maximum(np.nanstd(x,0),1e-6)
  e=d[target].to_numpy(float);e=(e-e.mean())/e.std();ys=d[PRIMARY_Y+SENS_Y].to_numpy(float);ys=(ys-ys.mean(0))/ys.std(0)
  take=d.CM_state.eq(state).to_numpy();groups=d.loc[take,'block200'].to_numpy()
  xx=prior.center(x[take],groups);ee=prior.center(e[take],groups).ravel();yy=prior.center(ys[take],groups)
  keep=xx.std(0)>.015;xx=xx[:,keep]
  er=ee-xx@np.linalg.lstsq(xx,ee,rcond=1e-8)[0];yr=yy-xx@np.linalg.lstsq(xx,yy,rcond=1e-8)[0]
  beta=er@yr/(er@er)
  got=coefficients[coefficients.stage.eq(stage)&coefficients.source_kind.eq(kind)&coefficients.source_layer.eq(layer)&coefficients.target.eq(target)&coefficients.CM_state.eq(state)&coefficients.geometry.eq('proximity50')&coefficients.adjustment.eq('within_region')].set_index('outcome')
  assert len(got)==len(beta)
  for k,v in zip(PRIMARY_Y+SENS_Y,beta):assert abs(v-got.loc[k,'coefficient'])<2e-6
  regressions+=len(beta)
 primary=coefficients.primary_analysis
 assert np.allclose(coefficients.loc[primary,'primary_grid_BH_q'],bh(coefficients.loc[primary,'within_section_p']),atol=1e-8)
 assert coefficients.loc[~primary,'primary_grid_BH_q'].isna().all()
 assert np.isfinite(coefficients[['coefficient','CI_low','CI_high']]).all().all()
 assert (coefficients.CI_low<=coefficients.CI_high).all()
 assert diagnostics.condition_number.max()<1e5
 val={'status':'PASS','independent_neighborhood_program_means_checked':checked,'independent_regression_coefficients_checked':regressions,
  'coefficient_rows':len(coefficients),'primary_hypotheses':int(primary.sum()),'global_primary_BH_recalculated':True,'missing_reference_classes_not_printed_as_zeros':True,
  'FB_program_genes_disjoint_from_CM_maturation_panel':True,'largest_design_condition_number':diagnostics.condition_number.max(),
  'same_CM_views_within_each_geometry':True,'new_biological_replicates':0,'segmentation_uncertainty_fully_propagated':False}
 (OUT/'numerical_validation.json').write_text(json.dumps(val,indent=2)+'\n');print('REFIT_VALIDATION_PASS',val,flush=True);export_inputs()

if __name__=='__main__':
 import argparse
 p=argparse.ArgumentParser();p.add_argument('--step',choices=['paired','spatial','validate','all'],default='all');a=p.parse_args()
 if a.step in ['paired','all']:paired()
 if a.step in ['spatial','all']:spatial()
 if a.step in ['validate','all']:validation()
