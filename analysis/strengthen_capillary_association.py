"""Predeclared P7 geometry, regional influence, subtype and local-match checks.

Tests challenge the existing association; no outcome-guided threshold selection,
new segmentation, imputed RNA, causal mechanism or animal replication is added.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','2')
os.environ.setdefault('MPLCONFIGDIR','/tmp/capillary-robustness-mpl')
from pathlib import Path
import json,hashlib
import numpy as np
import pandas as pd
from scipy import sparse,linalg
from scipy.spatial import cKDTree
from scipy.sparse.csgraph import dijkstra
from scipy.stats import t,false_discovery_control
from sklearn.preprocessing import SplineTransformer
import cardiac_subtype_niches as prior

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/capillary_association_robustness';OUT.mkdir(parents=True,exist_ok=True)
SEED=281028
def read(p,**kw):return pd.read_csv(ROOT/p,sep='\t',**kw)
def save(d,p,index=False):d.to_csv(OUT/p,sep='\t',index=index,float_format='%.9g',na_rep='NA')
def standard(x):return (x-x.mean(0))/np.maximum(x.std(0),1e-9)
def input_manifest():
 paths=['results/capillary_maturation_retest/P7_CM_retest_values.tsv.gz',
  'results/cardiac_spatial_pipeline/P7/atlas_annotations.tsv.gz',
  'results/cardiac_spatial_reconstruction/P7/vascular_interface_and_matrix_context.tsv.gz',
  'results/cardiac_spatial_reconstruction/P7/tissue_screened_Delaunay_graph.npz',
  'results/cardiac_spatial_pipeline/P7/proximity50_graph.npz']
 save(pd.DataFrame([dict(path=p,bytes=(ROOT/p).stat().st_size,
  sha256=hashlib.sha256((ROOT/p).read_bytes()).hexdigest()) for p in paths]),'input_checksums.tsv')

def model(d,exposure,outcomes,grid=200,offset=(0,0),extra=None,wild=False):
 """Within-region/type regression, clustered t and null-imposed wild scores."""
 cols=['capture','cell_area','distance_to_large_void_um','distance_to_outer_um',exposure]+outcomes
 if extra:cols+=extra
 valid=np.isfinite(d[cols].to_numpy(float)).all(1);z=d.loc[valid].copy()
 if len(z)<200:return [],None
 yy=standard(z[outcomes].to_numpy(float));ee=standard(z[[exposure]].to_numpy(float)).ravel()
 spline=SplineTransformer(n_knots=4,degree=3,include_bias=False).fit_transform(np.log1p(z[['capture','cell_area']]))
 cov=np.column_stack([spline,z.nucleus_missing.astype(float),np.log1p(z[['distance_to_large_void_um','distance_to_outer_um']]),pd.get_dummies(z.morphology_region,dtype=float)])
 if extra:cov=np.column_stack([cov,standard(z[extra].to_numpy(float))])
 groups=((z.x_centroid-offset[0])//grid).astype(int).astype(str)+'_'+((z.y_centroid-offset[1])//grid).astype(int).astype(str)+'::'+z.CM_state
 xx=prior.center(standard(cov),groups);ee=prior.center(ee,groups).ravel();yy=prior.center(yy,groups)
 _,r,piv=linalg.qr(xx,mode='economic',pivoting=True);rank=np.sum(abs(np.diag(r))>max(abs(r[0,0]),1e-12)*1e-9);xx=xx[:,piv[:rank]]
 er=ee-xx@np.linalg.lstsq(xx,ee,rcond=1e-10)[0]
 yr=yy-xx@np.linalg.lstsq(xx,yy,rcond=1e-10)[0]
 if er.std()<.01:return [],None
 den=er@er;beta=er@yr/den;resid=yr-er[:,None]*beta
 block=(z.x_centroid//400).astype(int).astype(str)+'_'+(z.y_centroid//400).astype(int).astype(str)
 codes=pd.factorize(block)[0];nb=codes.max()+1
 if nb<8:return [],None
 cluster_scores=np.zeros((nb,len(outcomes)));null_scores=np.zeros_like(cluster_scores);xxsum=np.zeros(nb);xysum=np.zeros_like(cluster_scores)
 for b in range(nb):
  take=codes==b;cluster_scores[b]=er[take]@resid[take];null_scores[b]=er[take]@yr[take]
  xxsum[b]=er[take]@er[take];xysum[b]=er[take]@yr[take]
 se=np.sqrt((cluster_scores**2).sum(0)/den**2*nb/(nb-1)*(len(z)-1)/(len(z)-rank-1))
 ci=t.ppf(.975,nb-1)*se
 wildp=np.full(len(outcomes),np.nan)
 if wild:
  rng=np.random.default_rng(SEED);B=9999;sign=rng.choice([-1.,1.],size=(B,nb))
  null=sign@null_scores/den
  wildp=(1+(abs(null)>=abs(beta)).sum(0))/(B+1)
 rows=[dict(exposure=exposure,outcome=o,grid_um=grid,offset_x=offset[0],offset_y=offset[1],
  CMs=len(z),spatial_blocks=nb,beta_SD=beta[j],CI_low=beta[j]-ci[j],CI_high=beta[j]+ci[j],
  within_section_p=2*t.sf(abs(beta[j]/max(se[j],1e-15)),nb-1),wild_cluster_null_p=wildp[j],
  exposure_residual_SD=er.std(),extra_controls=';'.join(extra or [])) for j,o in enumerate(outcomes)]
 # Hold fitted nuisance residuals fixed while removing each cluster: influence diagnostic, not refit validation.
 leave=(xysum.sum(0)[None,:]-xysum)/(den-xxsum[:,None])
 details=dict(index=z.index,codes=codes,er=er,yr=yr,beta=beta,leave=leave,groups=groups)
 return rows,details

def main():
 d=read('results/capillary_maturation_retest/P7_CM_retest_values.tsv.gz',index_col=0);d=d[d.strict].copy()
 f=read('results/cardiac_spatial_pipeline/P7/atlas_annotations.tsv.gz',index_col=0);assert f.index.is_unique
 b=read('results/cardiac_spatial_reconstruction/P7/vascular_interface_and_matrix_context.tsv.gz',index_col=0)
 d['near_boundary']=-np.log1p(b.loc[d.index,'QC_capillary__minimum_boundary_distance_um'])
 xy=f[['x_centroid','y_centroid']].to_numpy();ci=f.index.get_indexer(d.index)
 cap=f.atlas_identity.eq('Endothelial')&f.program_source_eligible&f.EC_nuclear_type.eq('Capillary-like')
 pubqc=f.atlas_identity.eq('Endothelial')&f.program_source_eligible&f.author_atlas_cell_state.eq('cap.EC')
 eligible=f.program_source_eligible
 masks={'cell_supported_capillary':pubqc,'Fibroblast':f.atlas_identity.eq('Fibroblast')&eligible,
  'noncapillary_EC':f.atlas_identity.eq('Endothelial')&eligible&~cap,
  'Pericyte':f.atlas_identity.eq('Pericyte')&eligible,
  'any_nonCM':f.atlas_identity.isin(['Endothelial','Fibroblast','Pericyte','Smooth muscle','Immune','Epicardial'])&eligible}
 for name,mask in masks.items():d['near_'+name]=-np.log1p(cKDTree(xy[mask]).query(xy[ci])[0])
 g=sparse.load_npz(ROOT/'results/cardiac_spatial_reconstruction/P7/tissue_screened_Delaunay_graph.npz').tocoo()
 lengths=np.linalg.norm(xy[g.row]-xy[g.col],axis=1)
 weighted=sparse.csr_matrix((lengths,(g.row,g.col)),shape=g.shape)
 path=dijkstra(weighted,directed=False,indices=np.flatnonzero(cap),min_only=True)[ci]
 d['path_um']=path;d['near_tissue_path']=-np.log1p(path)
 G=sparse.load_npz(ROOT/'results/cardiac_spatial_pipeline/P7/proximity50_graph.npz').tocsr()[ci]
 d['log_local_cells']=np.log1p(np.asarray(G@f.atlas_resolved.to_numpy(float)).ravel())
 d['log_capillary_count']=np.log1p(np.asarray(G@cap.to_numpy(float)).ravel())
 d['log_EC_capture']=np.asarray(G@(cap.to_numpy(float)*f.source_log_capture.fillna(0).to_numpy())).ravel()/np.maximum(np.asarray(G@cap.to_numpy(float)).ravel(),1)
 d['boundary_um']=b.loc[d.index,'QC_capillary__minimum_boundary_distance_um']
 base=['core_maturation21','native_maturation21','core_metabolic','core_sarcomere','core_excitation_contraction','core_cycling']
 rows=[];influence=[]
 primary=['near_QC_capillary','near_boundary','near_tissue_path']
 for exposure in primary+['near_cell_supported_capillary']:
  rr,info=model(d,exposure,base,wild=exposure in primary)
  rows += [dict(analysis='geometry',**r) for r in rr]
  if info:
   for j,o in enumerate(base):influence.append(dict(exposure=exposure,outcome=o,spatial_blocks=info['leave'].shape[0],
    minimum_leave_cluster_slope=info['leave'][:,j].min(),maximum_leave_cluster_slope=info['leave'][:,j].max(),
    positive_fraction=(info['leave'][:,j]>0).mean(),nuisance_residuals_held_fixed=True))
 for exposure in primary:
  for grid in [100,200,400]:
   for frac in [(0,0),(.5,0),(0,.5),(.5,.5)]:
    rr,_=model(d,exposure,['core_maturation21'],grid=grid,offset=(grid*frac[0],grid*frac[1]))
    rows += [dict(analysis='grid_sensitivity',**r) for r in rr]
  for state,z in d.groupby('CM_state'):
   rr,_=model(z,exposure,['core_maturation21','native_maturation21'])
   rows += [dict(analysis='CM_subtype',CM_subtype=state,**r) for r in rr]
 # Additional conditional questions are kept separate from the primary proximity effect.
 for extra in [['log_local_cells','log_EC_capture'],['log_local_cells','log_EC_capture','log_capillary_count'],
               ['log_local_cells','near_Fibroblast','near_noncapillary_EC','near_Pericyte']]:
  rr,_=model(d,'near_QC_capillary',base,extra=extra)
  rows += [dict(analysis='specificity_or_density',**r) for r in rr]
 # Same-depth cohort sensitivity retains a measured score after exactly 50 molecules.
 rr,_=model(d[d.core_depth50.notna()],'near_boundary',['core_maturation21','core_depth50'])
 rows += [dict(analysis='equal_depth_boundary',**r) for r in rr]
 result=pd.DataFrame(rows);primarymask=result.analysis.eq('geometry')&result.exposure.isin(primary)&result.outcome.eq('core_maturation21')
 result['primary_analysis']=primarymask;result['primary_BH_q']=np.nan
 result.loc[primarymask,'primary_BH_q']=false_discovery_control(result.loc[primarymask,'wild_cluster_null_p'])
 save(result,'P7_robustness_regressions.tsv');save(pd.DataFrame(influence),'P7_regional_influence.tsv')
 save(d[['CM_state','x_centroid','y_centroid','capture','cell_area','nucleus_missing','distance_QC_capillary','boundary_um','path_um','core_maturation21','native_maturation21','near_cell_supported_capillary','log_local_cells','log_capillary_count']], 'P7_added_geometry_values.tsv.gz',True)
 # Prespecified near/far thresholds, local disjoint pairing and no score-based match selection.
 x=np.column_stack([np.log1p(d.capture),np.log1p(d.cell_area)]);x=standard(x)
 high=d.distance_QC_capillary.le(25).to_numpy();low=d.distance_QC_capillary.ge(50).to_numpy()
 tree=cKDTree(d[['x_centroid','y_centroid']]);neighbors=tree.query_ball_point(d[['x_centroid','y_centroid']],100)
 rng=np.random.default_rng(SEED);hi=np.flatnonzero(high);rng.shuffle(hi);used=set();pairs=[]
 region=((d.x_centroid//200).astype(int).astype(str)+'_'+(d.y_centroid//200).astype(int).astype(str)).to_numpy()
 for i in hi:
  js=np.array([j for j in neighbors[i] if low[j] and j not in used and region[j]==region[i]
    and d.CM_state.iloc[j]==d.CM_state.iloc[i] and d.morphology_region.iloc[j]==d.morphology_region.iloc[i]],dtype=int)
  if len(js)==0:continue
  diff=np.abs(x[js]-x[i]);good=(diff<=.25).all(1);js=js[good];diff=diff[good]
  if len(js)==0:continue
  j=js[np.argmin(np.sum(diff**2,axis=1))];used.add(j)
  pairs.append(dict(near_CM=d.index[i],far_CM=d.index[j],CM_subtype=d.CM_state.iloc[i],
    near_distance_um=d.distance_QC_capillary.iloc[i],far_distance_um=d.distance_QC_capillary.iloc[j],
    physical_pair_distance_um=np.linalg.norm(d[['x_centroid','y_centroid']].iloc[i]-d[['x_centroid','y_centroid']].iloc[j]),
    log_capture_difference=x[i,0]-x[j,0],log_area_difference=x[i,1]-x[j,1]))
 pairs=pd.DataFrame(pairs);save(pairs,'P7_local_near_far_pairs.tsv.gz')
 assert len(pairs)>30 and pairs.near_CM.is_unique and pairs.far_CM.is_unique and set(pairs.near_CM).isdisjoint(pairs.far_CM)
 a=d.loc[pairs.near_CM];b=d.loc[pairs.far_CM];block=(a.x_centroid//400).astype(int).astype(str)+'_'+(a.y_centroid//400).astype(int).astype(str)
 codes=pd.factorize(block)[0];nb=codes.max()+1;boot=[];rng=np.random.default_rng(SEED)
 weights=np.array([np.bincount(rng.integers(nb,size=nb),minlength=nb) for _ in range(1999)])
 for outcome in base+['core_depth50']:
  delta=a[outcome].to_numpy()-b[outcome].to_numpy();valid=np.isfinite(delta)
  sums=np.bincount(codes[valid],weights=delta[valid],minlength=nb);counts=np.bincount(codes[valid],minlength=nb)
  sampled=(weights@sums)/np.maximum(weights@counts,1);sd=d[outcome].std(ddof=0)
  lo,hi=np.quantile(sampled/sd,[.025,.975])
  boot.append(dict(outcome=outcome,matched_pairs=int(valid.sum()),blocks=int(np.unique(codes[valid]).size),
    standardized_difference=delta[valid].mean()/sd,CI_low=lo,CI_high=hi,raw_difference=delta[valid].mean()))
 save(pd.DataFrame(boot),'P7_local_matched_outcomes.tsv')
 (OUT/'methods.json').write_text(json.dumps(dict(primary_tests=3,stage='P7',original_strict_CMs=len(d),
  tissue_path_reachable_CMs=int(np.isfinite(path).sum()),wild_cluster_draws=9999,cluster_um=400,
  spatial_null='Rademacher sign flips of null-imposed clustered partial-regression scores; not spatial relocation of vessels',
  matching='<=25 versus >=50um nearest capillary centroid; within100um, same200um grid, sameCMstate/morphology; logcapture and logarea differences <=.25 global SD; unique cells; no outcome matching',
  graph_path='Shortest sum of centroid edge lengths on frozen tissue-screened Delaunay graph; not true tissue geodesic or perfusion',
  boundary='Previous native CM boundary sampled at64points versus QC-supported capillary polygons; approximate2D distance',
  regional_influence='Leave one400um cluster out of fitted partial-regression sums; nuisance residuals fixed, not independent held-out validation',
  new_segmentation=False,new_biological_replicates=0,causality_claimed=False),indent=2)+'\n')
 input_manifest()
 print(result[result.primary_analysis].to_string(index=False));print(pd.DataFrame(boot).to_string(index=False))
if __name__=='__main__':main()
