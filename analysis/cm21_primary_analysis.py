"""The sole maturation outcome is the directly measured 21-gene Xenium index.

Maps retain every reviewed ventricular CM passing frozen QC. No anatomical
selection, published-label concordance, missing-gene imputation or additional
CM gene program defines the outcome. Existing non-CM identity labels are inputs.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','2')
os.environ.setdefault('MPLCONFIGDIR','/tmp/cm21-primary-mpl')
from pathlib import Path
import json,hashlib
import numpy as np
import pandas as pd
import h5py
from anndata.io import read_elem,sparse_dataset
from scipy import sparse,linalg
from scipy.spatial import cKDTree
from scipy.stats import norm,t,false_discovery_control
from sklearn.preprocessing import SplineTransformer
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/CM21_primary_spatial';OUT.mkdir(parents=True,exist_ok=True)
STAGES=['P0','P7','P14','P21'];SEED=281029
PANEL=pd.read_csv(ROOT/'results/fibroblast_reassessment/CM_maturation_Figure4b_panel.tsv',sep='\t').gene.tolist()
BROAD=['Endothelial','Fibroblast','Pericyte','Smooth muscle','Immune','Epicardial']
EC=['Capillary-like','Endocardial-like','Arterial-like','Venous-like','Lymphatic-like'];FB=['Postn','Gfpt2','Col8a1','Lsamp','Kirrel3']
INPUTS=set()
def read(p,**kw):INPUTS.add(str(p));return pd.read_csv(ROOT/p,sep='\t',**kw)
def save(d,p,index=False):d.to_csv(OUT/p,sep='\t',index=index,float_format='%.9g',na_rep='NA')
def z(x):return (x-np.mean(x,axis=0))/np.maximum(np.std(x,axis=0),1e-9)

def source_masks(stage,f):
 eligible=f.program_source_eligible.fillna(False);masks={b:f.atlas_identity.eq(b)&eligible for b in BROAD};catalog=[]
 for name in EC:masks['EC:'+name]=masks['Endothelial']&f.EC_nuclear_type.eq(name)
 fb=read(f'results/fibroblast_multimethod_identity/{stage}_multimethod_FB_annotations.tsv.gz',index_col=0)
 assess=read('results/cardiac_fibroblast_CM_refit/fibroblast_assignment_coverage.tsv')
 assess=assess[assess.stage.eq(stage)&assess.label_layer.eq('candidate')].set_index('family').reference_assessable
 for name in FB:
  mask=pd.Series(False,index=f.index);mask.loc[fb.index]=fb.program_source_eligible&fb.recommended_candidate.eq(name+'-associated')
  masks['FB:'+name]=mask&masks['Fibroblast']
 for name,m in masks.items():
  available=bool(assess.loc[name[3:]]) if name.startswith('FB:') else True
  catalog.append(dict(stage=stage,source=name,source_profiles=int(m.sum()),reference_assessable=available,
   support_status='Available' if available and m.sum()>=30 else 'Reference not assessable' if not available else 'Insufficient sources'))
 return masks,pd.DataFrame(catalog)

def pair_high_low(d,view='core'):
 """No source identity, distance or neighbor count participates in matching."""
 rank=d[view+'_stage_rank'];high=rank.ge(.75);low=rank.le(.25)
 cov=z(np.column_stack([np.log1p(d.CM_outcome_excluded_depth),np.log1p(d.cell_area)]))
 xy=d[['x_centroid','y_centroid']].to_numpy();grid=(xy//200).astype(int)
 groups=pd.Series(grid[:,0].astype(str)+'_'+grid[:,1].astype(str)+'::'+d.CM_state.to_numpy(),index=d.index)
 rows=[];rng=np.random.default_rng(SEED+int(d.stage.iloc[0][1:])+(0 if view=='core' else 100))
 for _,ids in groups.groupby(groups).groups.items():
  idx=d.index.get_indexer(ids);hi=idx[high.iloc[idx]];lo=idx[low.iloc[idx]]
  if not len(hi) or not len(lo):continue
  used=set()
  for i in rng.permutation(hi):
   candidates=np.array([j for j in lo if j not in used and d.nucleus_missing.iloc[j]==d.nucleus_missing.iloc[i]],dtype=int)
   if not len(candidates):continue
   diff=abs(cov[candidates]-cov[i]);physical=np.linalg.norm(xy[candidates]-xy[i],axis=1)
   keep=(diff.max(1)<=.25)&(physical<=100)
   if not keep.any():continue
   candidates=candidates[keep];diff=diff[keep];physical=physical[keep]
   k=np.argmin(np.square(diff).sum(1)+.05*np.square(physical/100));j=candidates[k];used.add(j)
   rows.append(dict(stage=d.stage.iloc[0],view=view,high_CM=d.index[i],low_CM=d.index[j],CM_state=d.CM_state.iloc[i],
    physical_separation_um=physical[k],log_capture_SD_difference=cov[i,0]-cov[j,0],log_area_SD_difference=cov[i,1]-cov[j,1],
    block400=f'{int(xy[i,0]//400)}_{int(xy[i,1]//400)}',
    high_score=d[view+'_score21'].iloc[i],low_score=d[view+'_score21'].iloc[j]))
 return pd.DataFrame(rows)

def paired_effect(values,pairs,d,metric):
 hi=d.index.get_indexer(pairs.high_CM);lo=d.index.get_indexer(pairs.low_CM)
 # Positive distance contrast means higher-score CMs are closer; presence is high minus low.
 delta=values[lo]-values[hi] if metric=='closer_um' else values[hi]-values[lo]
 codes=pd.factorize(pairs.block400)[0];nb=codes.max()+1
 occupied=np.bincount(codes,weights=((values[hi]>0)|(values[lo]>0)).astype(float),minlength=nb)
 if metric=='presence_pp' and ((occupied>=3).sum()<8 or occupied.sum()<50):return None
 sums=np.bincount(codes,weights=delta,minlength=nb);counts=np.bincount(codes,minlength=nb)
 rng=np.random.default_rng(SEED);weights=np.array([np.bincount(rng.integers(nb,size=nb),minlength=nb) for _ in range(1999)])
 boot=(weights@sums)/np.maximum(weights@counts,1);point=delta.mean();loCI,hiCI=np.quantile(boot,[.025,.975]);se=boot.std(ddof=1)
 factor=100 if metric=='presence_pp' else 1
 return dict(metric=metric,matched_pairs=len(pairs),spatial_blocks=nb,estimate=factor*point,CI_low=factor*loCI,CI_high=factor*hiCI,
  within_section_p=2*norm.sf(abs(point)/max(se,1e-12)),high_mean=float(values[hi].mean())*(factor if metric=='presence_pp' else 1),
  low_mean=float(values[lo].mean())*(factor if metric=='presence_pp' else 1))

def regression(d,exposure,view='core',local=False,anatomy=False):
 """Supporting continuous-score check, retaining unclassified anatomy."""
 y=z(d[[view+'_score21']].to_numpy()).ravel();e=z(exposure)
 C=SplineTransformer(n_knots=4,degree=3,include_bias=False).fit_transform(np.log1p(d[['CM_outcome_excluded_depth','cell_area']]))
 C=np.column_stack([C,d.nucleus_missing.astype(float),pd.get_dummies(d.CM_state,dtype=float)])
 if anatomy:C=np.column_stack([C,pd.get_dummies(d.anatomy_region,dtype=float)])
 C=z(C)
 if local:
  groups=(d.x_centroid//200).astype(int).astype(str)+'_'+(d.y_centroid//200).astype(int).astype(str)+'::'+d.CM_state
  def center(x):
   df=pd.DataFrame(np.atleast_2d(x).T if x.ndim==1 else x,index=d.index);return (df-df.groupby(groups).transform('mean')).to_numpy()
  C=center(C);e=center(e).ravel();y=center(y).ravel()
 else:C=np.column_stack([np.ones(len(d)),C])
 _,r,piv=linalg.qr(C,mode='economic',pivoting=True);rank=np.sum(abs(np.diag(r))>max(abs(r[0,0]),1e-12)*1e-9);C=C[:,piv[:rank]]
 er=e-C@np.linalg.lstsq(C,e,rcond=1e-10)[0]
 if er.std()<.01:return None
 yr=y-C@np.linalg.lstsq(C,y,rcond=1e-10)[0];den=er@er;beta=er@yr/den;res=yr-er*beta
 block=(d.x_centroid//400).astype(int).astype(str)+'_'+(d.y_centroid//400).astype(int).astype(str)
 score=pd.Series(er*res,index=d.index).groupby(block).sum().to_numpy();ng=len(score)
 se=np.sqrt(np.square(score).sum()/den**2*ng/(ng-1)*(len(d)-1)/(len(d)-rank-1));ci=t.ppf(.975,ng-1)*se
 return dict(beta_SD=beta,CI_low=beta-ci,CI_high=beta+ci,within_section_p=2*t.sf(abs(beta/max(se,1e-12)),ng-1),CMs=len(d),blocks=ng)

def main():
 assert len(PANEL)==len(set(PANEL))==21
 cal=read('results/cardiac_spatial_reconstruction/score_reference_calibration.tsv').set_index(['view','scope','gene']).sort_index()
 cohort=[];detection=[];descriptive=[];paired=[];models=[];allpairs=[];coverage=[];genechecks=[]
 for stage in STAGES:
  print('CM21_PRIMARY',stage,flush=True)
  f=read(f'results/cardiac_spatial_pipeline/{stage}/atlas_annotations.tsv.gz',index_col=0)
  eligible=f.ventricular_CM_strict&f.atlas_identity.eq('Cardiomyocyte')&f.pass_expression_QC&f.core_transcripts.ge(20)&f.CM_outcome_excluded_depth.ge(20)
  d=f.loc[eligible,['stage','x_centroid','y_centroid','cell_area','nucleus_area','CM_outcome_excluded_depth','author_atlas_cell_state']].copy()
  d['nucleus_missing']=d.nucleus_area.isna();known=['v.CM','Ankrd1.CM','Slit2.CM','Ddc.CM','p.CM']
  d['CM_state']=d.author_atlas_cell_state.where(d.author_atlas_cell_state.isin(known),'CM subtype unresolved')
  score_old=read(f'results/CM_maturation21_nonCM_context/{stage}/CM_maturation21_scores.tsv.gz',index_col=0)
  assert d.index.isin(score_old.index).all()
  for view in ['core','native']:
   path=ROOT/(f'results/xenium_nuclear_anchored_segmentation/{stage}/nuclear_anchored_counts.h5ad' if view=='core' else f'results/xenium_ventricular_mapping/{stage}/{stage}_xenium_counts_and_candidate_mapping.h5ad')
   INPUTS.add(str(path))
   with h5py.File(path) as h:
    ids=read_elem(h['obs']).index;genes=read_elem(h['var']).index;ii=ids.get_indexer(d.index);jj=genes.get_indexer(PANEL)
    assert (ii>=0).all() and (jj>=0).all() and len(genes)==403
    X=sparse_dataset(h['layers/high_confidence'] if view=='core' else h['X'])[ii].tocsr();counts=X[:,jj].toarray().astype(float)
    total=np.asarray(X.sum(1)).ravel();del X
   ref=cal.loc[(view,'all')].reindex(PANEL);log=np.log1p(10000*counts/np.maximum(total,1)[:,None]);zz=(log-ref.reference_mean.to_numpy())/ref.reference_SD.to_numpy()
   score=zz.mean(1);error=np.max(abs(score-score_old.loc[d.index,view+'_Z__maturation21']))
   assert error<2e-6
   genechecks.append(dict(stage=stage,view=view,genes=21,CMs=len(d),max_score_reconstruction_error=error))
   d[view+'_score21']=score;d[view+'_stage_rank']=pd.Series(score,index=d.index).rank(pct=True)
   d[view+'_score_group']=np.where(d[view+'_stage_rank'].le(.25),'Lower quartile',np.where(d[view+'_stage_rank'].ge(.75),'Higher quartile','Middle half'))
   d[view+'_detected_panel_genes']=(counts>0).sum(1);d[view+'_panel_molecules']=counts.sum(1)
   for j,g in enumerate(PANEL):detection.append(dict(stage=stage,view=view,gene=g,CMs=len(d),fraction_detected=(counts[:,j]>0).mean(),mean_log1p_CP10k=log[:,j].mean()))
  # Anatomical labels are only a covariate sensitivity; no CM is removed by them.
  annotations=read(f'results/morphology_stratified_nonCM/{stage}_CM_region_annotations_and_scores.tsv.gz',index_col=0)
  d['anatomy_region']=annotations.region.reindex(d.index).fillna('Unclassified')
  # The expanded cohort adds previously excluded cells; their anatomy remains unclassified rather than guessed.
  cohort.append(dict(stage=stage,QC_eligible_ventricular_CMs=len(d),old_restricted_CMs=int(score_old.strict_ventricular_map_cohort.sum()),
   added_CMs=len(d)-int(score_old.strict_ventricular_map_cohort.sum()),lower_quartile=int(d.core_score_group.eq('Lower quartile').sum()),
   higher_quartile=int(d.core_score_group.eq('Higher quartile').sum()),missing_nuclear_area=int(d.nucleus_missing.sum()),
   median_measured_21genes=int(d.core_detected_panel_genes.median())))
  masks,catalog=source_masks(stage,f);xy=f[['x_centroid','y_centroid']].to_numpy();ci=f.index.get_indexer(d.index)
  G=sparse.load_npz(ROOT/f'results/cardiac_spatial_pipeline/{stage}/proximity50_graph.npz').tocsr()[ci]
  INPUTS.add(f'results/cardiac_spatial_pipeline/{stage}/proximity50_graph.npz')
  pairsets={view:pair_high_low(d,view) for view in ['core','native']}
  for view,pairs in pairsets.items():
   assert len(pairs)>100 and pairs.high_CM.is_unique and pairs.low_CM.is_unique and set(pairs.high_CM).isdisjoint(pairs.low_CM)
   allpairs.append(pairs)
  for name,mask in masks.items():
   row=catalog[catalog.source.eq(name)].iloc[0].to_dict()
   assert not d.index.isin(f.index[mask]).any(),(stage,name,'source/receiver conflict')
   if row['support_status']!='Available':coverage.append(row);continue
   distance=cKDTree(xy[mask]).query(xy[ci])[0];count=np.asarray(G@mask.to_numpy(float)).ravel();presence=count>0
   near=distance<=100;nearblocks=np.unique((d.loc[near,['x_centroid','y_centroid']].to_numpy()//400).astype(int),axis=0)
   row.update(CMs_within100um=int(near.sum()),nearby_blocks=len(nearblocks),primary_supported=near.sum()>=100 and len(nearblocks)>=8)
   coverage.append(row)
   for view in ['core','native']:
    scoregroup=d[view+'_score_group']
    for label in ['Lower quartile','Higher quartile']:
     take=scoregroup.eq(label).to_numpy()
     descriptive.append(dict(stage=stage,view=view,source=name,score_group=label,CMs=int(take.sum()),
      mean_nearest_um=distance[take].mean(),median_nearest_um=np.median(distance[take]),fraction_with_source_50um=presence[take].mean()))
    if not row['primary_supported']:continue
    pairs=pairsets[view]
    for metric,values in [('closer_um',distance),('presence_pp',presence.astype(float))]:
     effect=paired_effect(values,pairs,d,metric)
     if effect:paired.append(dict(stage=stage,view=view,source=name,**effect))
    for adjustment,local,anatomy in [('Capture_and_CM_identity',False,False),('Within_local_region',True,False),('Anatomy_sensitivity',True,True)]:
     result=regression(d,-np.log1p(distance),view,local,anatomy)
     if result:models.append(dict(stage=stage,view=view,source=name,adjustment=adjustment,**result))
   d['distance__'+name]=distance;d['presence50__'+name]=presence;d['neighbor_count50__'+name]=count
  d.index.name='cell_id';save(d,f'{stage}_CM21_scores_and_neighborhoods.tsv.gz',True)
  pairs=pairsets['core'];d.loc[d.index.isin(pairs.high_CM),'matched_group']='Higher';d.loc[d.index.isin(pairs.low_CM),'matched_group']='Lower'
 effect=pd.DataFrame(paired);primary=effect.view.eq('core')&effect.metric.eq('closer_um')
 effect['distance_comparison']=primary;effect['matched_distance_BH_q']=np.nan;effect.loc[primary,'matched_distance_BH_q']=false_discovery_control(effect.loc[primary,'within_section_p'])
 secondary=effect.view.eq('core')&effect.metric.eq('presence_pp');effect['neighborhood_presence_BH_q']=np.nan
 effect.loc[secondary,'neighborhood_presence_BH_q']=false_discovery_control(effect.loc[secondary,'within_section_p'])
 save(effect,'high_low_CM21_matched_comparisons.tsv')
 model=pd.DataFrame(models);sel=model.view.eq('core')&model.adjustment.eq('Within_local_region');model['primary_analysis']=sel;model['primary_BH_q']=np.nan
 model.loc[sel,'primary_BH_q']=false_discovery_control(model.loc[sel,'within_section_p']);save(model,'continuous_CM21_associations.tsv')
 save(pd.concat(allpairs,ignore_index=True),'high_low_CM21_matched_pairs.tsv.gz');save(pd.DataFrame(cohort),'CM21_stage_cohorts.tsv')
 save(pd.DataFrame(detection),'measured_21gene_detection.tsv');save(pd.DataFrame(descriptive),'all_CM21_high_low_neighborhood_descriptions.tsv')
 save(pd.DataFrame(coverage),'nonCM_source_coverage.tsv');save(pd.DataFrame(genechecks),'independent_21gene_reconstruction.tsv')
 save(pd.DataFrame(dict(gene=PANEL)),'CM_maturation_21gene_panel.tsv')
 (OUT/'methods.json').write_text(json.dumps(dict(primary_CM_outcome='One directly measured21-gene Xenium score; no other CM gene programs or projected genes',
  genes=PANEL,scoring='Mean of21 frozen gene-Z values after log1p CP10k normalization over403 measured targets; separate frozen core/native calibrations',
  exact_author_AddModuleScore_recovered=False,CM_QC='Reviewed ventricular CM; expressionQC pass; >=20core transcripts; >=20outcome-excluded counts; nucleus-free profiles retained',
  anatomical_exclusion=False,published_label_concordance_required=False,quartiles='Within each stage, across all QC-eligible ventricular CMs; rank ties averaged',
  primary_comparison='Continuous directly measured21-gene score across all QC-eligible CMs; nearness adjusted for capture, size, CM identity and local200um regions; no anatomical exclusion',
  matched_sensitivity='High versus low stage-quartile CMs matched within200um grid and available CM subtype; <=100um separation; capture/area calipers <=.25SD; same missing-nucleus flag; unique cells',
  matched_role_reason='Strict high/low matching retains few cells because score and capture are correlated; it is a sensitivity rather than a replacement for the all-CM question. This is exploratory reanalysis, not a prospectively registered trial.',
  matching_does_not_use_neighbors=True,positive_primary_effect='Higher21-gene CM score associated with source nearness; continuous coefficient in SD',
  positive_matched_distance_effect='Mean low-score distance minus high-score distance, micrometers; positive means higher-score CMs are closer',
  primary_tests=int(sel.sum()),matched_distance_tests=int(primary.sum()),presence_secondary_tests=int(secondary.sum()),spatial_bootstrap='1999 resamples of400um pair-anchor clusters; within-section uncertainty',
  source_IDs='Existing reviewed source labels; no new reference clustering or gene imputation',new_biological_replicates=0,only21_CM_outcome_genes=True),indent=2)+'\n')
 (OUT/'input_paths.json').write_text(json.dumps(sorted(INPUTS),indent=2)+'\n')
 print(model[model.primary_analysis].sort_values('primary_BH_q').head(20).to_string(index=False),flush=True)
if __name__=='__main__':main()
