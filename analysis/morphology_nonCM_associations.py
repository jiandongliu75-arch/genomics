"""Measured CM maturation versus non-CM subtypes within frozen image ROIs."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','2')
os.environ.setdefault('MPLCONFIGDIR','/tmp/anatomy-mpl')
from pathlib import Path
import json,hashlib
import numpy as np
import pandas as pd
from scipy import sparse,linalg
from scipy.sparse.csgraph import dijkstra
from scipy.spatial import cKDTree
from scipy.stats import t,false_discovery_control
from sklearn.preprocessing import SplineTransformer
from matplotlib.path import Path as PolygonPath
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/morphology_stratified_nonCM'
STAGES=['P0','P7','P14','P21'];BROAD=['Endothelial','Fibroblast','Pericyte','Smooth muscle','Immune','Epicardial']
EC=['Capillary-like','Endocardial-like','Arterial-like','Venous-like','Lymphatic-like'];FB=['Postn','Gfpt2','Col8a1','Lsamp','Kirrel3']
COMP=['maturation21','metabolic','sarcomere','excitation_contraction','cycling'];INPUTS=set()
def read(p,**kw):INPUTS.add(str(p));return pd.read_csv(ROOT/p,sep='\t',**kw)
def save(d,p,index=False):d.to_csv(OUT/p,sep='\t',index=index,float_format='%.9g',na_rep='NA')
def zscore(x):
 x=np.asarray(x,float);good=np.isfinite(x);out=np.full(len(x),np.nan)
 if good.any():out[good]=(x[good]-x[good].mean())/max(x[good].std(),1e-9)
 return out
def assign(stage,xy):
 obj=json.loads((OUT/'morphology_ROIs.geojson').read_text());a=pd.DataFrame(index=np.arange(len(xy)))
 a['ROI']='Unclassified';a['region']='Unclassified';a['territory']='Unresolved';a['confidence']='none'
 # Non-overlapping priority: compact, cavity-facing, then low-confidence candidates.
 for feature in obj['features']:
  p=feature['properties']
  if p['stage']!=stage:continue
  poly=np.array(feature['geometry']['coordinates'][0]);take=PolygonPath(poly,closed=True).contains_points(xy,radius=1e-8)
  a.loc[take,['ROI','region','territory','confidence']]=[p['ROI_id'],p['region'],p['image_territory'],p['confidence']]
 return a
def fit(d,target,outcomes,scope,context=True,extra=None):
 fields=[target,'capture','cell_area','distance_to_large_void_um','distance_to_outer_um']+outcomes+(extra or [])
 take=np.isfinite(d[fields].to_numpy(float)).all(1);z=d.loc[take]
 if len(z)<200:return [],'Fewer than200 finite CMs'
 blocks=(z.x_centroid//400).astype(int).astype(str)+'_'+(z.y_centroid//400).astype(int).astype(str)
 code=pd.factorize(blocks)[0];nb=code.max()+1
 if nb<8:return [],'Fewer than8 spatial clusters'
 if 'source_local_support' in z:
  supported=z.source_local_support.to_numpy(bool)
  occupied=np.bincount(code[supported],minlength=nb)
  if supported.sum()<100 or (occupied>=5).sum()<8:return [],'Insufficient local source support across8 clusters'
 C=SplineTransformer(n_knots=4,degree=3,include_bias=False).fit_transform(np.log1p(z[['capture','cell_area']]))
 C=np.column_stack([C,z.nucleus_missing.astype(float),np.log1p(z[['distance_to_large_void_um','distance_to_outer_um']])])
 if extra:C=np.column_stack([C,z[extra].to_numpy(float)])
 C=(C-C.mean(0))/np.maximum(C.std(0),1e-8)
 group=(z.x_centroid//200).astype(int).astype(str)+'_'+(z.y_centroid//200).astype(int).astype(str)+'::'+z.CM_state
 if context:group=group+'::'+z.ROI
 def center(x):
  x=np.asarray(x);x=x[:,None] if x.ndim==1 else x
  inverse=pd.factorize(group)[0];counts=np.bincount(inverse);means=np.zeros((len(counts),x.shape[1]));np.add.at(means,inverse,x)
  return x-means[inverse]/counts[inverse,None]
 C=center(C);e=center(z[target].to_numpy()).ravel();Y=center(z[outcomes].to_numpy())
 _,r,piv=linalg.qr(C,mode='economic',pivoting=True);rank=np.sum(abs(np.diag(r))>max(abs(r[0,0]),1e-12)*1e-9);C=C[:,piv[:rank]]
 er=e-C@np.linalg.lstsq(C,e,rcond=1e-10)[0]
 if er.std()<.015:return [],'Insufficient residual exposure variation'
 yr=Y-C@np.linalg.lstsq(C,Y,rcond=1e-10)[0];den=er@er;beta=er@yr/den;res=yr-er[:,None]*beta
 scores=np.zeros((nb,len(outcomes)))
 for j in range(nb):scores[j]=er[code==j]@res[code==j]
 se=np.sqrt((scores*scores).sum(0)/den**2*nb/(nb-1)*(len(z)-1)/(len(z)-rank-1));ci=t.ppf(.975,nb-1)*se
 return [dict(scope=scope,exposure=target,outcome=o,beta_SD=beta[j],CI_low=beta[j]-ci[j],CI_high=beta[j]+ci[j],
  within_section_p=2*t.sf(abs(beta[j]/max(se[j],1e-15)),nb-1),CMs=len(z),spatial_blocks=nb,
  residual_exposure_SD=er.std(),context_ROI_conditioned=context,extra_controls=';'.join(extra or [])) for j,o in enumerate(outcomes)],'Fitted'

def sources(stage,f):
 masks={};meta=[];eligible=f.program_source_eligible.fillna(False)
 for b in BROAD:masks['Broad:'+b]=f.atlas_identity.eq(b)&eligible
 for name in EC:masks['EC:'+name]=masks['Broad:Endothelial']&f.EC_nuclear_type.eq(name)
 fb=read(f'results/fibroblast_multimethod_identity/{stage}_multimethod_FB_annotations.tsv.gz',index_col=0)
 ref=read('results/cardiac_fibroblast_CM_refit/fibroblast_assignment_coverage.tsv');ref=ref[ref.stage.eq(stage)&ref.label_layer.eq('candidate')].set_index('family')
 for name in FB:
  m=pd.Series(False,index=f.index);m.loc[fb.index]=fb.program_source_eligible&fb.recommended_candidate.eq(name+'-associated')
  masks['FB:'+name]=m
 masks['EC:Cycling-associated']=masks['Broad:Endothelial']&f.author_atlas_cell_state.eq('p.EC')&f.EC_cycling.gt(0)
 # Alternative EC source rules are sensitivity layers, not replacements.
 for name,public in [('Capillary-like','cap.EC'),('Endocardial-like','endocardial.EC'),('Arterial-like','arteriole.EC'),('Lymphatic-like','lym.EC')]:
  masks['EC_public:'+name]=masks['Broad:Endothelial']&f.author_atlas_cell_state.eq(public)
 for key,m in masks.items():
  assessable=bool(ref.loc[key[3:],'reference_assessable']) if key.startswith('FB:') else True
  meta.append(dict(stage=stage,source=key,eligible_source_profiles=int(m.sum()),reference_assessable=assessable,
   source_layer='Public label plus reviewed broad EC' if key.startswith('EC_public:') else 'Current reviewed source definitions',
   subtype_status='RNA-supported candidate' if key.startswith(('FB:','EC:')) else 'Broad lineage' if key.startswith('Broad:') else 'Alternative candidate'))
 return masks,pd.DataFrame(meta)

def main():
 freeze=json.loads((OUT/'anatomy_freeze.json').read_text());assert hashlib.sha256((OUT/'morphology_ROIs.geojson').read_bytes()).hexdigest()==freeze['mask_sha256']
 results=[];coverage=[];cohorts=[];skips=[];descriptive=[]
 for stage in STAGES:
  print('ANATOMY_NONCM',stage,flush=True)
  d=read(f'results/capillary_maturation_retest/{stage}_CM_retest_values.tsv.gz',index_col=0);d=d[d.strict].copy()
  f=read(f'results/cardiac_spatial_pipeline/{stage}/atlas_annotations.tsv.gz',index_col=0)
  a=assign(stage,d[['x_centroid','y_centroid']].to_numpy());a.index=d.index;d=d.join(a)
  d['annotated_primary']=d.confidence.isin(['high','moderate'])
  # Frozen stage scales make region-specific standardized slopes comparable.
  outcomes=[]
  for view in ['core','native']:
   for comp in COMP:
    name=view+'_'+comp;d['Y__'+name]=zscore(d[name]);outcomes.append('Y__'+name)
  for group,z in d.groupby('region'):
   cohorts.append(dict(stage=stage,region=group,CMs=len(z),ROIs=z.ROI.nunique(),
    spatial_blocks_400=z.assign(block=(z.x_centroid//400).astype(int).astype(str)+'_'+(z.y_centroid//400).astype(int).astype(str)).block.nunique(),
    primary_eligible=group in ['Compact wall','Cavity-facing myocardium']))
   for comp in COMP:descriptive.append(dict(stage=stage,region=group,outcome=comp,CMs=len(z),mean_core_score=z['core_'+comp].mean()))
  masks,catalog=sources(stage,f);xy=f[['x_centroid','y_centroid']].to_numpy();ci=f.index.get_indexer(d.index)
  graph_path=ROOT/f'results/cardiac_spatial_reconstruction/{stage}/tissue_screened_Delaunay_graph.npz'
  INPUTS.add(str(graph_path.relative_to(ROOT)));g=sparse.load_npz(graph_path).tocoo()
  weighted=sparse.csr_matrix((np.linalg.norm(xy[g.row]-xy[g.col],axis=1),(g.row,g.col)),shape=g.shape)
  G=sparse.load_npz(ROOT/f'results/cardiac_spatial_pipeline/{stage}/proximity50_graph.npz').tocsr()[ci]
  INPUTS.add(f'results/cardiac_spatial_pipeline/{stage}/proximity50_graph.npz')
  d['density']=zscore(np.log1p(np.asarray(G@f.atlas_resolved.to_numpy(float)).ravel()))
  for key,m in masks.items():
   row=catalog[catalog.source.eq(key)].iloc[0].to_dict();n=int(m.sum())
   if not row['reference_assessable'] or n<30:
    coverage.append(dict(**row,analysis_status='Reference not assessable' if not row['reference_assessable'] else 'Fewer than30 sources'));continue
   distance=cKDTree(xy[m]).query(xy[ci])[0]
   path=dijkstra(weighted,directed=False,indices=np.flatnonzero(m),min_only=True)[ci]
   near=(distance<=100)&d.annotated_primary.to_numpy();nb=len(np.unique((d.loc[near,['x_centroid','y_centroid']].to_numpy()//400).astype(int),axis=0))
   row.update(annotated_CMs_within100um=int(near.sum()),nearby_spatial_blocks=nb,reachable_CMs=int(np.isfinite(path).sum()))
   if near.sum()<100 or nb<8:
    coverage.append(dict(**row,analysis_status='Insufficient local CM coverage'));continue
   coverage.append(dict(**row,analysis_status='Tested'))
   exposures={'path':zscore(-np.log1p(path)),'centroid':zscore(-np.log1p(distance))}
   count=np.asarray(G@m.to_numpy(float)).ravel();d['source_capture']=zscore(np.asarray(G@(m.to_numpy(float)*f.source_log_capture.fillna(0).to_numpy(float))).ravel()/np.maximum(count,1))
   denom=(masks['Broad:Fibroblast'] if key.startswith('FB:') else masks['Broad:Endothelial'] if key.startswith(('EC:','EC_public:')) else f.atlas_resolved)
   fraction=count/np.maximum(np.asarray(G@denom.to_numpy(float)).ravel(),1)
   exposures['fraction50']=zscore(fraction)
   for metric,value in exposures.items():
    d['E']=value
    d['source_local_support']=((path<=100) if metric=='path' else (distance<=100) if metric=='centroid' else count>0)
    # Primary comparison is among annotated image regions; all other scopes are explicit sensitivities.
    scopes={'Annotated_regions':d.annotated_primary,'Compact_wall':d.region.eq('Compact wall'),
     'Cavity_facing':d.region.eq('Cavity-facing myocardium')}
    if metric=='path':scopes['All_frozen_CMs']=pd.Series(True,index=d.index)
    for scope,take in scopes.items():
     rr,status=fit(d.loc[take],'E',outcomes,scope,context=scope!='All_frozen_CMs')
     results += [dict(stage=stage,source=key,metric=metric,**r) for r in rr]
     if not rr:skips.append(dict(stage=stage,source=key,metric=metric,scope=scope,reason=status))
    if metric=='path':
     rr,_=fit(d.loc[d.annotated_primary],'E',outcomes,'Density_capture_sensitivity',extra=['density','source_capture'])
     results += [dict(stage=stage,source=key,metric=metric,**r) for r in rr]
     # Named CM-state interactions remain secondary and require sufficient spatial support.
     for state,take in d[d.annotated_primary].groupby('CM_state'):
      if key not in ['EC:Capillary-like','EC:Endocardial-like','FB:Postn','FB:Gfpt2','FB:Col8a1']:continue
      rr,status=fit(take,'E',outcomes,'CM_state:'+state)
      results += [dict(stage=stage,source=key,metric=metric,**r) for r in rr]
   d=d.drop(columns=['E'])
  # Measured fibroblast programs complement uncertain family labels at late ages.
  pr=read('results/cardiac_fibroblast_CM_refit/fibroblast_per_profile_programs.tsv.gz')
  pr=pr[pr.stage.eq(stage)&pr.view.eq('core')].set_index('cell_id');good=pr.informative
  si=f.index.get_indexer(pr.index);gb=G[:,si];den=np.asarray(gb@good.to_numpy(float)).ravel()
  d['FB_program_capture']=zscore(np.asarray(gb@(np.log1p(pr.view_depth.to_numpy())*good.to_numpy())).ravel()/np.maximum(den,1))
  d['FB_program_count']=zscore(np.log1p(den))
  for program in [c for c in pr if c.startswith('program__')]:
   d['E']=zscore(np.asarray(gb@(pr[program].to_numpy()*good.to_numpy())).ravel()/np.maximum(den,1))
   d['source_local_support']=den>0
   rr,_=fit(d.loc[d.annotated_primary&(den>0)],'E',outcomes,'Measured_FB_programs_without_source_controls')
   results += [dict(stage=stage,source=program,metric='mean_program50',**r) for r in rr]
   rr,_=fit(d.loc[d.annotated_primary&(den>0)],'E',outcomes,'Measured_FB_programs',extra=['FB_program_capture','FB_program_count'])
   results += [dict(stage=stage,source=program,metric='mean_program50',**r) for r in rr]
  save(d.drop(columns=['E','source_capture','source_local_support'],errors='ignore'),f'{stage}_CM_region_annotations_and_scores.tsv.gz',True)
  save(catalog,f'{stage}_source_catalog.tsv')
 result=pd.DataFrame(results)
 primary=(result.scope.eq('Annotated_regions')&result.metric.eq('path')&~result.source.str.startswith('EC_public:')&result.outcome.eq('Y__core_maturation21'))
 result['primary_analysis']=primary;result['primary_BH_q']=np.nan
 result.loc[primary,'primary_BH_q']=false_discovery_control(result.loc[primary,'within_section_p'])
 # Components, subtypes and alternative scopes are secondary; all are retained and corrected together.
 secondary=~primary&result.outcome.str.startswith('Y__core_');result['secondary_BH_q']=np.nan
 result.loc[secondary,'secondary_BH_q']=false_discovery_control(result.loc[secondary,'within_section_p'])
 save(result,'nonCM_maturation_associations.tsv');save(pd.DataFrame(coverage),'source_coverage_and_assessability.tsv')
 save(pd.DataFrame(cohorts),'anatomical_CM_coverage.tsv');save(pd.DataFrame(descriptive),'region_score_descriptions.tsv');save(pd.DataFrame(skips),'unassessable_region_tests.tsv')
 (OUT/'association_design.json').write_text(json.dumps(dict(anatomy_sha256=freeze['mask_sha256'],primary_tests=int(primary.sum()),
  primary='Within annotated compact/cavity-facing patches; tissue-screened weighted graph nearness; overall21 score; core RNA; all ages and supported sources corrected together',
  comparators='Centroid nearness and tissue-screened50um neighbor fractions; fixed stage scales; nativeRNA; source density/capture; CM subtypes and separate regions',
  broad_only_nonCM=['Immune','Pericyte','Smooth muscle','Epicardial'],fine_immune_taxonomy_available=False,
  fibroblast_programs='Previously measured rank-matched programs, not imputed RNA or fate labels',
  spatial_uncertainty='400um cluster CR1 t intervals; conditioning within ROI x CM subtype x200um grid',
  anatomical_labels_are_candidates=True,LV_RV_verified=False,papillary_verified=False,large_vessels_verified=False,
  new_segmentation=False,new_biological_replicates=0,causality=False),indent=2)+'\n')
 (OUT/'analysis_input_paths.json').write_text(json.dumps(sorted(INPUTS),indent=2)+'\n')
 print(result[result.primary_analysis].sort_values('primary_BH_q').head(20).to_string(index=False),flush=True)
if __name__=='__main__':main()
