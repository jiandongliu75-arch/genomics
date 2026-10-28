"""Independent count-space, geometry and explicit-regression verification."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','2')
from pathlib import Path
import json,hashlib
import numpy as np
import pandas as pd
import h5py
from anndata.io import read_elem,sparse_dataset
from scipy import sparse
from scipy.stats import t
from sklearn.preprocessing import SplineTransformer
from shapely import Polygon,STRtree,line_interpolate_point

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/CM21_paperstyle_extension'
def read(p,**kw):return pd.read_csv(p,sep='\t',**kw)
def bh(p):
 p=np.asarray(p);order=np.argsort(p);q=np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]
 ans=np.empty(len(p));ans[order]=np.minimum(q,1);return ans
def main():
 panel=read(OUT/'measured_21gene_panel.tsv').gene.tolist();assert len(set(panel))==21
 weights=read(OUT/'score_gene_weights.tsv.gz');models=read(OUT/'spatial_models.tsv');inter=read(OUT/'fibroblast_context_interactions.tsv')
 cal=read(ROOT/'results/cardiac_spatial_reconstruction/score_reference_calibration.tsv').set_index(['view','scope','gene']).sort_index()
 audit=[];geocheck=[];modelcheck=[]
 for stage in ['P0','P7','P14','P21']:
  print('VALIDATE',stage,flush=True)
  d=read(OUT/(stage+'_CM21_extension_values.tsv.gz'),index_col=0)
  base=read(ROOT/f'results/CM21_primary_spatial/{stage}_CM21_scores_and_neighborhoods.tsv.gz',index_col=0)
  assert d.index.equals(base.index)
  f=read(ROOT/f'results/cardiac_spatial_pipeline/{stage}/atlas_annotations.tsv.gz',index_col=0)
  cap=(f.atlas_identity.eq('Endothelial')&f.program_source_eligible&f.EC_nuclear_type.eq('Capillary-like'))
  G=sparse.load_npz(ROOT/f'results/cardiac_spatial_pipeline/{stage}/proximity50_graph.npz').tocsr()
  assert np.array_equal(np.asarray(G[f.index.get_indexer(d.index)]@cap.to_numpy(float)).ravel(),d.capillary_count50)
  rng=np.random.default_rng(771);sample=rng.choice(len(d),size=384,replace=False)
  for view in ['native','core']:
   p=ROOT/(f'results/xenium_ventricular_mapping/{stage}/{stage}_xenium_counts_and_candidate_mapping.h5ad' if view=='native' else f'results/xenium_nuclear_anchored_segmentation/{stage}/nuclear_anchored_counts.h5ad')
   with h5py.File(p) as h:
    ids=read_elem(h['obs']).index;genes=read_elem(h['var']).index
    x=sparse_dataset(h['X'] if view=='native' else h['layers/high_confidence'])[ids.get_indexer(d.index[sample])].toarray().astype(float)
   total=x.sum(1);log=np.log1p(10000*x/np.maximum(total,1)[:,None])
   w=weights[weights.stage.eq(stage)&weights.view.eq(view)].pivot(index='gene',columns='score',values='weight').reindex(genes)
   assert np.max(abs(w.sum(0)))<2e-7
   error=np.max(abs(log@w.module_ensemble.to_numpy()-d.iloc[sample][view+'_module21'].to_numpy()))
   assert error<3e-7,(stage,view,error)
   assert np.max(abs(total-x[:,genes.get_indexer(panel)].sum(1)-d.iloc[sample][view+'_excluded21_depth'].to_numpy()))<.01
   means=cal.loc[(view,'all')].reindex(genes).reference_mean.to_numpy();bins=np.empty(len(genes),int)
   bins[np.argsort(means,kind='stable')]=np.arange(len(genes))*8//len(genes)
   ww=[]
   for k in range(50):
    rng2=np.random.default_rng(20261028+k);selected=[]
    for g in panel:selected.extend(rng2.choice(np.flatnonzero(bins==bins[genes.get_loc(g)]),20,replace=False))
    selected=np.unique(selected);v=np.zeros(403);v[genes.get_indexer(panel)]=1/21;v[selected]-=1/len(selected);ww.append(v)
    assert np.max(abs(v-w['draw_'+str(k)].to_numpy()))<1e-8
   assert np.max(abs(np.mean(ww,axis=0)-w.module_ensemble.to_numpy()))<1e-8
   for g in panel:assert abs(w.loc[g,'omit_'+g])<1e-15,(g,'omitted gene retained')
   audit.append(dict(stage=stage,view=view,random_count_profiles=384,max_independent_score_error=error,all50_control_draws_reconstructed=True,all21_omissions_remove_gene_from_both_module_and_controls=True))
  # Explicit partial least squares with independently constructed fixed-effect removal.
  n=len(d);y=d.native_module21.to_numpy();y=(y-y.mean())/y.std()
  e=-np.log1p(base['distance__EC:Capillary-like'].to_numpy());e=(e-e.mean())/e.std()
  C=SplineTransformer(n_knots=4,degree=3,include_bias=False).fit_transform(np.log1p(np.column_stack([d.native_excluded21_depth,base.cell_area])))
  C=np.column_stack([C,base.nucleus_missing.astype(float),pd.get_dummies(base.CM_state,dtype=float)])
  C=(C-C.mean(0))/np.maximum(C.std(0),1e-12)
  groups=(base.x_centroid//200).astype(int).astype(str)+'_'+(base.y_centroid//200).astype(int).astype(str)+'::'+base.CM_state
  codes=pd.factorize(groups)[0];ng=codes.max()+1
  def demean(a):
   df=pd.DataFrame(a);return (df-df.groupby(codes).transform('mean')).to_numpy()
  C=demean(C);e=demean(e)[:,0];y=demean(y)[:,0]
  er=e-C@np.linalg.lstsq(C,e,rcond=1e-9)[0];yr=y-C@np.linalg.lstsq(C,y,rcond=1e-9)[0]
  beta=er@yr/(er@er);res=yr-beta*er
  blocks=pd.factorize((base.x_centroid//400).astype(int).astype(str)+'_'+(base.y_centroid//400).astype(int).astype(str))[0];nb=blocks.max()+1
  sums=np.bincount(blocks,weights=er*res);rank=np.linalg.matrix_rank(C,tol=1e-8)
  se=np.sqrt(np.square(sums).sum()/(er@er)**2*nb/(nb-1)*(n-1)/(n-rank-ng-1));ci=t.ppf(.975,nb-1)*se
  actual=models[models.stage.eq(stage)&models.vascular_primary&models.exposure.eq('capillary_nearness')].iloc[0]
  assert abs(beta-actual.beta_SD)<2e-7
  assert max(abs(beta-ci-actual.CI_low),abs(beta+ci-actual.CI_high))<2e-7
  modelcheck.append(dict(stage=stage,independent_native_capillary_beta=beta,reported_beta=actual.beta_SD,absorbed_region_state_groups=ng,blocks=nb,independent_CI_low=beta-ci,independent_CI_high=beta+ci))
  # Re-read original boundaries for 24 CMs and all eligible capillary polygons.
  pos=np.flatnonzero(d.capillary_apposition2.to_numpy()>0);neg=np.flatnonzero(d.capillary_apposition2.to_numpy()==0)
  chosen=np.r_[rng.choice(pos,min(12,len(pos)),replace=False),rng.choice(neg,min(12,len(neg)),replace=False)]
  wanted=set(f.index[cap])|set(d.index[chosen]);parts=[]
  p=next(Path('/workspace/heart-atlas/data/xenium').glob('WT_'+('P07' if stage=='P7' else stage)+'__*'))/'cell_boundaries.csv.gz'
  for chunk in pd.read_csv(p,usecols=['cell_id','vertex_x','vertex_y'],chunksize=400000):
   parts.append(chunk[chunk.cell_id.isin(wanted)])
  v=pd.concat(parts,ignore_index=True);polys={cid:Polygon(q[['vertex_x','vertex_y']].to_numpy()) for cid,q in v.groupby('cell_id',sort=False) if len(q)>=3}
  tree=STRtree([polys[c] for c in f.index[cap] if c in polys and not polys[c].is_empty])
  for cid in d.index[chosen]:
   points=line_interpolate_point(polys[cid].boundary,(np.arange(64)+.5)/64,normalized=True)
   _,dist=tree.query_nearest(points,return_distance=True,all_matches=False)
   error=abs((dist<=2).mean()-d.loc[cid,'capillary_apposition2']);assert error<1e-8
   assert abs(dist.min()-d.loc[cid,'capillary_boundary_um'])<1e-5
   geocheck.append(dict(stage=stage,cell_id=cid,independent_fraction_2um=(dist<=2).mean(),reported_fraction_2um=d.loc[cid,'capillary_apposition2']))
 take=models.vascular_primary;assert take.sum()==12;assert np.max(abs(bh(models.loc[take,'wild_score_p'])-models.loc[take,'vascular_BH_q']))<1e-8
 # Eight-significant-digit TSV serialization can shift recomputed q by ~1e-8.
 take=models.nonCM_BH_q.notna();assert take.sum()==48;assert np.max(abs(bh(models.loc[take,'within_section_p'])-models.loc[take,'nonCM_BH_q']))<3e-8
 take=inter.view.eq('native')&inter.score.eq('module_ensemble');assert take.sum()==12;assert np.max(abs(bh(inter.loc[take,'within_section_p'])-inter.loc[take,'context_BH_q']))<1e-8
 label=read(OUT/'source_label_sensitivity.tsv');basechecks=0
 for r in label.itertuples():
  assert r.sources>=30 and np.isfinite([r.beta_SD,r.CI_low,r.CI_high]).all()
  if r.metric!='nearness' or r.layer not in ['reviewed_nuclear_supported','recommended_candidate']:continue
  match=models[models.stage.eq(r.stage)&models.view.eq(r.view)&models.source.eq(r.source)&models.score.eq('module_ensemble')&models.exposure.eq('source_nearness')]
  if len(match):
   assert len(match)==1 and abs(match.iloc[0].beta_SD-r.beta_SD)<2e-7
   assert abs(match.iloc[0].CI_low-r.CI_low)<2e-7
   basechecks+=1
 for name,rows in [('independent_count_validation.tsv',audit),('independent_boundary_validation.tsv',geocheck),('independent_model_validation.tsv',modelcheck)]:pd.DataFrame(rows).to_csv(OUT/name,sep='\t',index=False)
 receipt=dict(status='PASS',all_four_cohorts_preserved=True,count_views_checked=8,random_count_profiles_checked=3072,
  independently_rebuilt_control_draws=400,one_gene_omission_weight_checks=168,original_boundary_profiles_checked=len(geocheck),
  explicit_local_regressions_checked=4,BH_families_verified=[12,48,12],source_label_baselines_checked=basechecks,new_biological_replicates=0)
 (OUT/'validation.json').write_text(json.dumps(receipt,indent=2)+'\n');print(receipt,flush=True)
if __name__=='__main__':main()
