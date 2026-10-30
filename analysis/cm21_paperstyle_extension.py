"""Exploratory, measured-21-gene-only score and spatial-context challenges.

Frozen settings are written before reading new expression outcomes. This is a
paper-style reconstruction: original control genes/masks/graph are unavailable.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','2')
os.environ.setdefault('OMP_NUM_THREADS','2')
from pathlib import Path
import json, hashlib
import numpy as np
import pandas as pd
import h5py
from anndata.io import read_elem, sparse_dataset
from scipy import sparse, linalg
from scipy.stats import t, spearmanr, false_discovery_control
from scipy.sparse.csgraph import dijkstra
from sklearn.preprocessing import SplineTransformer
import cm21_primary_analysis as old

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/CM21_paperstyle_extension'
STAGES=old.STAGES; PANEL=old.PANEL; SEED=20261028
CONFIG=dict(genes=PANEL,cohort='Unchanged all-QC ventricular CM cohort; no anatomical exclusion',
 primary_score='Mean of 50 AddModuleScore-style draws: 8 expression bins, 20 controls per module gene, all 403 measured genes eligible, unique union per draw',
 score_reference='Frozen pooled all-profile log1p CP10k means; fixed expression-rank bins across stages within each segmentation view',
 score_implementation='Rank-balanced expression bins (not exact Seurat cut_number); equal module-gene weights minus mean of unique control union; controls may include module genes, as in default Seurat pool',
 main_view='Deposited whole-cell Xenium counts; confidence-supported core counts are a segmentation sensitivity, not nuclear-only truth',
 settings=[dict(name='primary_8bins_20controls',bins=8,controls=20,exclude_panel=False),
 dict(name='literal100_3bins',bins=3,controls=100,exclude_panel=False),
 dict(name='fine24_5controls',bins=24,controls=5,exclude_panel=False),
 dict(name='exclude21_8bins_20controls',bins=8,controls=20,exclude_panel=True)],
 draws=50,seed=SEED,primary_radius_um=50,boundary_apposition_um=2,
 primary_vascular_metrics=['capillary_nearness','capillary_count50','capillary_apposition2'],
 primary_uncertainty='400um clustered t; null-imposed Rademacher clustered-score test, 9999 draws, for vascular tests; not vessel relocation',
 primary_vascular_correction='BH across all supported stage x metric tests, native primary score, within-local model; wild-score p',
 nonCM_correction='Separate BH family: all supported native primary-score subtype/broad source nearness tests',
 context_correction='Separate BH family: native capillary-nearness x fibroblast count50 interaction tests',
 adjustments='View-specific non-panel molecule count and cell area cubic splines; missing nuclear area; available CM state; 200um grid x CM-state fixed effects',
 context_model='Standardized capillary nearness, log1p FB-family count50, their product; local total resolved neighbor count included',
 anatomy='Existing candidate morphology labels sensitivity only; unknown anatomy retained; no exact author masks',
 score_sensitivity='Every one-gene omission, all 50 control draws, alternative bin/control pools, unchanged earlier gene-Z score',
 limitations=['One section per stage; no independent-heart replication','No exact author score, graph or mask recovery',
 'No perfusion, oxygen, fate or causal effect inferred','No additional CM outcome genes or imputation','Existing native boundaries; no new segmentation claimed'])

def read(p,**kw): return pd.read_csv(ROOT/p,sep='\t',**kw)
def save(d,p,index=False): d.to_csv(OUT/p,sep='\t',index=index,float_format='%.8g',na_rep='NA')
def z(a):
 a=np.asarray(a,float);return (a-a.mean(0))/np.maximum(a.std(0),1e-12)
def center(a,codes):
 a=np.asarray(a,float);one=a.ndim==1
 if one:a=a[:,None]
 n=np.bincount(codes); sums=np.zeros((len(n),a.shape[1]));np.add.at(sums,codes,a)
 b=a-sums[codes]/n[codes,None]
 return b[:,0] if one else b

class Design:
 def __init__(self,d,view,grid=200,offset=0,local=True,anatomy=False,cluster=400):
  self.d=d
  a=SplineTransformer(n_knots=4,degree=3,include_bias=False).fit_transform(np.log1p(d[[view+'_excluded21_depth','cell_area']]))
  a=np.column_stack([a,d.nucleus_missing.astype(float),pd.get_dummies(d.CM_state,dtype=float)])
  if anatomy:a=np.column_stack([a,pd.get_dummies(d.anatomy_region,dtype=float)])
  self.groups=None
  self.absorbed=0
  if local:
   g=((d.x_centroid-offset)//grid).astype(int).astype(str)+'_'+((d.y_centroid-offset)//grid).astype(int).astype(str)+'::'+d.CM_state
   self.groups=pd.factorize(g)[0]
   self.absorbed=int(self.groups.max()+1)
  a=z(a)
  if self.groups is not None:a=center(a,self.groups)
  else:a=np.column_stack([np.ones(len(d)),a])
  q,r,p=linalg.qr(a,mode='economic',pivoting=True)
  self.rank=int((abs(np.diag(r))>max(abs(r[0,0]),1e-12)*1e-9).sum());self.Q=q[:,:self.rank]
  b=(d.x_centroid//cluster).astype(int).astype(str)+'_'+(d.y_centroid//cluster).astype(int).astype(str)
  self.blocks=pd.factorize(b)[0];self.nb=self.blocks.max()+1
 def residual(self,a,standard=True):
  a=z(a) if standard else np.asarray(a,float)
  if self.groups is not None:a=center(a,self.groups)
  return a-self.Q@(self.Q.T@a)
 def fit(self,x,Y,names,extra=None,wild=False):
  er=self.residual(x);yr=self.residual(Y)
  if yr.ndim==1:yr=yr[:,None]
  rank=self.rank
  if extra is not None:
   ee=self.residual(extra);ee=np.atleast_2d(ee).T if ee.ndim==1 else ee
   q,r,p=linalg.qr(ee,mode='economic',pivoting=True)
   k=int((abs(np.diag(r))>max(abs(r[0,0]),1e-12)*1e-9).sum());q=q[:,:k]
   er-=q@(q.T@er);yr-=q@(q.T@yr);rank+=k
  den=er@er
  if er.std()<.005:return []
  beta=er@yr/den;res=yr-er[:,None]*beta
  cs=np.zeros((self.nb,yr.shape[1]));np.add.at(cs,self.blocks,er[:,None]*res)
  se=np.sqrt((cs**2).sum(0)/den**2*self.nb/(self.nb-1)*(len(er)-1)/(len(er)-rank-self.absorbed-1))
  ci=t.ppf(.975,self.nb-1)*se;p=2*t.sf(abs(beta/np.maximum(se,1e-15)),self.nb-1)
  wp=np.full(len(beta),np.nan)
  if wild:
   ns=np.zeros_like(cs);np.add.at(ns,self.blocks,er[:,None]*yr)
   rng=np.random.default_rng(SEED);hits=np.zeros(len(beta),int)
   for start in range(0,9999,333):
    n=min(333,9999-start);null=rng.choice([-1.,1.],size=(n,self.nb))@ns/den
    hits+=(abs(null)>=abs(beta)).sum(0)
   wp=(1+hits)/10000
  return [dict(score=names[j],CMs=len(er),blocks=self.nb,beta_SD=beta[j],CI_low=beta[j]-ci[j],CI_high=beta[j]+ci[j],
   within_section_p=p[j],wild_score_p=wp[j],exposure_residual_SD=er.std()) for j in range(len(beta))]

def control_weights(genes,means,setting,draw,omit=None):
 rng=np.random.default_rng(SEED+draw)
 allowed=np.array([g not in PANEL for g in genes]) if setting['exclude_panel'] else np.ones(len(genes),bool)
 if omit is not None:allowed &= genes.to_numpy()!=omit
 # Binning is fixed before exclusions/omissions, preserving matched expression strata.
 order=np.argsort(means,kind='stable');bins=np.empty(len(genes),int)
 bins[order]=np.minimum(np.arange(len(genes))*setting['bins']//len(genes),setting['bins']-1)
 chosen=[];module=[g for g in PANEL if g!=omit]
 for g in module:
  j=genes.get_loc(g);candidates=np.flatnonzero((bins==bins[j])&allowed)
  if len(candidates)<setting['controls']:raise ValueError((setting,g,len(candidates)))
  chosen.extend(rng.choice(candidates,size=setting['controls'],replace=False))
 union=np.unique(chosen);w=np.zeros(len(genes));w[genes.get_indexer(module)]=1/len(module);w[union]-=1/len(union)
 return w,union,bins

def main():
 OUT.mkdir(parents=True,exist_ok=True);(OUT/'analysis_settings.json').write_text(json.dumps(CONFIG,indent=2)+'\n')
 cal=read('results/cardiac_spatial_reconstruction/score_reference_calibration.tsv').set_index(['view','scope','gene'])
 weights=[];controlaudit=[];stability=[];models=[];context=[];geometry=[];rings=[];coverage=[];scoreq=[];manifest=[]
 for stage in STAGES:
  print('PAPERSTYLE',stage,flush=True)
  d=read(f'results/CM21_primary_spatial/{stage}_CM21_scores_and_neighborhoods.tsv.gz',index_col=0)
  f=read(f'results/cardiac_spatial_pipeline/{stage}/atlas_annotations.tsv.gz',index_col=0)
  masks,catalog=old.source_masks(stage,f);ci=f.index.get_indexer(d.index)
  G=sparse.load_npz(ROOT/f'results/cardiac_spatial_pipeline/{stage}/proximity50_graph.npz').tocsr()[ci]
  d['resolved_count50']=np.asarray(G@f.atlas_resolved.to_numpy(float)).ravel()
  cap=masks['EC:Capillary-like'].to_numpy();capcount=np.asarray(G@cap.astype(float)).ravel()
  b=read(f'results/cardiac_spatial_reconstruction/{stage}/vascular_interface_and_matrix_context.tsv.gz',index_col=0).reindex(d.index)
  d['capillary_apposition2']=b['QC_capillary__apposition_fraction_2um']
  d['capillary_boundary_um']=b['QC_capillary__minimum_boundary_distance_um']
  d['capillary_count50']=capcount;d['capillary_nearness']=-np.log1p(d['distance__EC:Capillary-like'])
  assert np.isfinite(d[['capillary_apposition2','capillary_boundary_um']]).all().all()
  # Physical path and graph-ring measures use the same existing tissue-screened graph.
  D=sparse.load_npz(ROOT/f'results/cardiac_spatial_reconstruction/{stage}/tissue_screened_Delaunay_graph.npz').tocoo()
  xy=f[['x_centroid','y_centroid']].to_numpy()
  W=sparse.csr_matrix((np.linalg.norm(xy[D.row]-xy[D.col],axis=1),(D.row,D.col)),shape=D.shape)
  path=dijkstra(W,directed=False,indices=np.flatnonzero(cap),min_only=True)[ci]
  ring=dijkstra(D.tocsr(),directed=False,unweighted=True,indices=np.flatnonzero(cap),min_only=True)[ci]
  d['capillary_tissue_path_um']=np.where(np.isfinite(path),path,np.nan)
  d['capillary_graph_ring']=np.where(np.isfinite(ring),np.minimum(ring,10),np.nan)
  coverage.append(dict(stage=stage,CMs=len(d),finite_path=int(np.isfinite(path).sum()),boundary_valid=int(d.capillary_apposition2.notna().sum()),capillary_sources=int(cap.sum())))
  counts_by_source={name:np.asarray(G@m.to_numpy(float)).ravel() for name,m in masks.items() if catalog.set_index('source').loc[name,'support_status']=='Available'}
  oldcoverage=read('results/CM21_primary_spatial/nonCM_source_coverage.tsv');supported=set(oldcoverage[oldcoverage.stage.eq(stage)&oldcoverage.primary_supported.fillna(False)].source)
  for view in ['native','core']:
   print(stage,view,'scores',flush=True)
   p=ROOT/(f'results/xenium_ventricular_mapping/{stage}/{stage}_xenium_counts_and_candidate_mapping.h5ad' if view=='native' else f'results/xenium_nuclear_anchored_segmentation/{stage}/nuclear_anchored_counts.h5ad')
   with h5py.File(p) as h:
    ids=read_elem(h['obs']).index;genes=read_elem(h['var']).index
    x=sparse_dataset(h['X'] if view=='native' else h['layers/high_confidence'])[ids.get_indexer(d.index)].tocsr()
    total=np.asarray(x.sum(1)).ravel();panelcounts=np.asarray(x[:,genes.get_indexer(PANEL)].sum(1)).ravel()
    log=np.log1p(1e4*x.toarray().astype(float)/np.maximum(total,1)[:,None]);del x
   d[view+'_excluded21_depth']=total-panelcounts
   means=cal.loc[(view,'all')].reindex(genes).reference_mean.to_numpy()
   # All control draws and omission scores are retained as small weight vectors.
   setting=CONFIG['settings'][0]
   ws=np.column_stack([control_weights(genes,means,setting,k)[0] for k in range(50)])
   seedscore=log@ws;primary=seedscore.mean(1)
   Y=[primary,d[view+'_score21'].to_numpy()];names=['module_ensemble','earlier_geneZ']
   for alternative in CONFIG['settings'][1:]:
    aw=np.column_stack([control_weights(genes,means,alternative,k)[0] for k in range(50)]).mean(1)
    Y.append(log@aw);names.append(alternative['name'])
    for j,g in enumerate(genes):weights.append(dict(stage=stage,view=view,score=alternative['name'],gene=g,weight=aw[j]))
   for omit in PANEL:
    ow=np.column_stack([control_weights(genes,means,setting,k,omit)[0] for k in range(50)]).mean(1)
    Y.append(log@ow);names.append('omit_'+omit)
    for j,g in enumerate(genes):weights.append(dict(stage=stage,view=view,score='omit_'+omit,gene=g,weight=ow[j]))
   for k in range(50):Y.append(seedscore[:,k]);names.append('draw_'+str(k))
   Y=np.column_stack(Y)
   for j,g in enumerate(genes):weights.append(dict(stage=stage,view=view,score='module_ensemble',gene=g,weight=ws[j].mean()))
   for k in range(50):
    _,union,bins=control_weights(genes,means,setting,k)
    controlaudit.append(dict(stage=stage,view=view,draw=k,unique_controls=len(union),module_genes_in_controls=len(set(genes[union])&set(PANEL)),minimum_bin_size=int(np.bincount(bins).min())))
    for j,g in enumerate(genes):weights.append(dict(stage=stage,view=view,score='draw_'+str(k),gene=g,weight=ws[j,k]))
   primaryrank=pd.Series(primary).rank(pct=True).to_numpy();high=primaryrank>=.75;low=primaryrank<=.25
   for j,name in enumerate(names):
    r=pd.Series(Y[:,j]).rank(pct=True).to_numpy();hh=r>=.75;ll=r<=.25
    stability.append(dict(stage=stage,view=view,score=name,Spearman_to_primary=spearmanr(primary,Y[:,j]).statistic,
     high_quartile_Jaccard=(hh&high).sum()/max((hh|high).sum(),1),low_quartile_Jaccard=(ll&low).sum()/max((ll|low).sum(),1)))
   d[view+'_module21']=primary;d[view+'_module21_rank']=primaryrank
   scoreq.append(dict(stage=stage,view=view,CMs=len(d),median=np.median(primary),p10=np.quantile(primary,.1),p90=np.quantile(primary,.9),
    raw_score_capture_Spearman=spearmanr(primary,total-panelcounts).statistic,core_native='Not a biological replicate'))
   designs={a:Design(d,view,local=a!='Across_section',anatomy=a=='Anatomy_sensitivity') for a in ['Across_section','Within_local','Anatomy_sensitivity']}
   for exposure in CONFIG['primary_vascular_metrics']:
    xx=d[exposure].to_numpy();xx=np.log1p(xx) if exposure=='capillary_count50' else xx
    for adjustment,design in designs.items():
     ys=Y if adjustment=='Within_local' else Y[:,:1];nn=names if adjustment=='Within_local' else names[:1]
     rr=design.fit(xx,ys,nn,wild=view=='native' and adjustment=='Within_local')
     models.extend(dict(stage=stage,view=view,source='EC:Capillary-like',exposure=exposure,adjustment=adjustment,**v) for v in rr)
    geometry.append(dict(stage=stage,view=view,exposure=exposure,raw_Spearman=spearmanr(primary,xx).statistic,
     low_quartile_mean=float(d.loc[low,exposure].mean()),high_quartile_mean=float(d.loc[high,exposure].mean()),
     low_quartile_nearest_capillary_um=d.loc[low,'distance__EC:Capillary-like'].mean(),
     high_quartile_nearest_capillary_um=d.loc[high,'distance__EC:Capillary-like'].mean()))
   for label,x in [('capillary_boundary_nearness',-np.log1p(d.capillary_boundary_um)),
                   ('capillary_tissue_nearness',-np.log1p(d.capillary_tissue_path_um)),
                   ('capillary_graph_nearness',-d.capillary_graph_ring)]:
    valid=np.isfinite(x.to_numpy());design=Design(d.loc[valid],view)
    rr=design.fit(x.to_numpy()[valid],Y[valid,:1],names[:1])
    models.extend(dict(stage=stage,view=view,source='EC:Capillary-like',exposure=label,adjustment='Alternative_geometry',**v) for v in rr)
   G30=sparse.load_npz(ROOT/f'results/cardiac_spatial_pipeline/{stage}/proximity30_graph.npz').tocsr()[ci]
   count30=np.asarray(G30@cap.astype(float)).ravel()
   rr=designs['Within_local'].fit(np.log1p(count30),Y[:,:1],names[:1])
   models.extend(dict(stage=stage,view=view,source='EC:Capillary-like',exposure='capillary_count30',adjustment='Radius_sensitivity',**v) for v in rr)
   # Density-specific question adds total resolved neighbors; all-cell coverage is a contextual quantity.
   for exposure in ['capillary_count50','capillary_nearness']:
    x=d[exposure].to_numpy();x=np.log1p(x) if exposure=='capillary_count50' else x
    rr=designs['Within_local'].fit(x,Y[:,:1],names[:1],extra=np.log1p(d.resolved_count50.to_numpy()))
    models.extend(dict(stage=stage,view=view,source='EC:Capillary-like',exposure=exposure,adjustment='Local_total_cell_count',**v) for v in rr)
   # Grid and larger uncertainty-block checks, without outcome-driven choice of scale.
   for grid,offset,cluster in [(100,0,400),(200,100,400),(400,0,400),(200,0,800)]:
    design=Design(d,view,grid=grid,offset=offset,cluster=cluster)
    for exposure in CONFIG['primary_vascular_metrics']:
     x=d[exposure].to_numpy();x=np.log1p(x) if exposure=='capillary_count50' else x
     rr=design.fit(x,Y[:,:1],names[:1]);models.extend(dict(stage=stage,view=view,source='EC:Capillary-like',exposure=exposure,adjustment=f'Grid{grid}_offset{offset}_cluster{cluster}',**v) for v in rr)
   # All available neighboring identities, without inventing late unavailable families.
   for source in sorted(supported):
    col='distance__'+source
    if col not in d:continue
    rr=designs['Within_local'].fit(-np.log1p(d[col]),Y[:,:5],names[:5])
    models.extend(dict(stage=stage,view=view,source=source,exposure='source_nearness',adjustment='Within_local',**v) for v in rr)
   for source in [s for s in sorted(supported) if s.startswith('FB:')]:
    count=counts_by_source[source];fb=z(np.log1p(count));near=z(d.capillary_nearness.to_numpy());interaction=near*fb
    rr=designs['Within_local'].fit(interaction,Y[:,:5],names[:5],extra=np.column_stack([near,fb,np.log1p(d.resolved_count50)]))
    context.extend(dict(stage=stage,view=view,FB_family=source,FB_sources=int(masks[source].sum()),CMs_with_FB50=int((count>0).sum()),**v) for v in rr)
    d['count50__'+source]=count
   for k in sorted(np.unique(d.capillary_graph_ring.dropna())):
    a=d.capillary_graph_ring.eq(k).to_numpy();rings.append(dict(stage=stage,view=view,ring=int(k),CMs=int(a.sum()),mean_module21=primary[a].mean(),median_module21=np.median(primary[a])))
   del log,Y,seedscore
  keep=['stage','x_centroid','y_centroid','CM_state','native_module21','native_module21_rank','core_module21','core_module21_rank','native_excluded21_depth','core_excluded21_depth','capillary_nearness','capillary_count50','resolved_count50','capillary_apposition2','capillary_boundary_um','capillary_tissue_path_um','capillary_graph_ring']+[c for c in d if c.startswith('count50__')]
  d.index.name='cell_id';save(d[keep],stage+'_CM21_extension_values.tsv.gz',True)
 m=pd.DataFrame(models);m['vascular_primary']=m.view.eq('native')&m.score.eq('module_ensemble')&m.adjustment.eq('Within_local')&m.exposure.isin(CONFIG['primary_vascular_metrics'])
 m['vascular_BH_q']=np.nan;m.loc[m.vascular_primary,'vascular_BH_q']=false_discovery_control(m.loc[m.vascular_primary,'wild_score_p'])
 take=m.view.eq('native')&m.score.eq('module_ensemble')&m.adjustment.eq('Within_local')&m.exposure.eq('source_nearness')
 m['nonCM_BH_q']=np.nan;m.loc[take,'nonCM_BH_q']=false_discovery_control(m.loc[take,'within_section_p'])
 c=pd.DataFrame(context);c['context_BH_q']=np.nan;take=c.view.eq('native')&c.score.eq('module_ensemble');c.loc[take,'context_BH_q']=false_discovery_control(c.loc[take,'within_section_p'])
 save(m,'spatial_models.tsv');save(c,'fibroblast_context_interactions.tsv');save(pd.DataFrame(stability),'score_stability.tsv')
 save(pd.DataFrame(weights),'score_gene_weights.tsv.gz');save(pd.DataFrame(controlaudit),'control_draw_audit.tsv');save(pd.DataFrame(geometry),'raw_vascular_comparisons.tsv')
 save(pd.DataFrame(rings),'tissue_graph_ring_means.tsv');save(pd.DataFrame(coverage),'geometry_coverage.tsv');save(pd.DataFrame(scoreq),'stage_score_summaries.tsv')
 save(pd.DataFrame(dict(gene=PANEL)),'measured_21gene_panel.tsv')
 print(m[m.vascular_primary].to_string(index=False),flush=True);print(c[c.view.eq('native')].to_string(index=False),flush=True)

if __name__=='__main__':main()
