"""Refit CM relationships with current FB candidates and directly measured programs.

Frozen coordinates/CM cohort/graphs; no new segmentation, classifier training,
spatial imputation, functional measurements or biological replicates.
"""
from pathlib import Path
import os,json,argparse,time,hashlib
os.environ.setdefault('OPENBLAS_NUM_THREADS','2')
import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse,linalg
from scipy.stats import norm
import cardiac_subtype_niches as prior
from fibroblast_multimethod_identity import refs,queries,INPUTS as FB_INPUTS
from fibroblast_multimethod_diagnostics import score_modules
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/cardiac_fibroblast_CM_refit'
STAGES=['P0','P7','P14','P21'];FAMILIES=['Postn','Gfpt2','Col8a1','Lsamp','Kirrel3']
PANELS={
 'Col8a1.Fib':['Cmklr1','Col8a1','Col6a3','Htra3'],
 'Gfpt2.Fib':['Gfpt2','Pi16','Adgrd1'],
 'Lsamp.Fib':['Col8a1','Itgbl1','Lsamp'],
 'Postn.Fib':['Col24a1','Ptn','Tnc','Postn'],
 'Fibrillar_matrix':['Col1a1','Col1a2','Col3a1'],
 'Basement_matrix':['Lama4','Lama5','Lamb1','Lamb2','Lamc1','Col4a1','Col4a2','Hspg2'],
 'Elastic_matrix':['Eln','Fbln1','Fbln2']}
MODULES=list(PANELS);GEOMETRIES=['proximity50','proximity30','contact']
LABELS={'candidate':'recommended_candidate','crosschecked':'crosschecked_candidate','legacy':'fine_identity'}
PRIMARY_Y=['core_Z__maturation21','core_Z__metabolic','core_Z__sarcomere','core_Z__excitation_contraction','core_Z__cycling']
SENS_Y=['core_Z__maturation20_noMyl2','core_Z__sarcomere_noMyl2']+['native_Z__'+k for k in ['maturation21','metabolic','sarcomere','excitation_contraction','cycling']]
SEED=261027;DRAWS=512;INPUTS=set()
def read(p,**kw):
 p=Path(p);p=p if p.is_absolute() else ROOT/p;INPUTS.add(str(p));return pd.read_csv(p,sep='\t',**kw)
def save(d,name):
 OUT.mkdir(parents=True,exist_ok=True);d.to_csv(OUT/name,sep='\t',index=False,na_rep='NA',float_format='%.9g')
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
def export_inputs():
 p=OUT/'input_paths.json';old=set(json.loads(p.read_text())) if p.exists() else set()
 extra={str(ROOT/k) if not Path(k).is_absolute() else k for k in FB_INPUTS}
 p.write_text(json.dumps(sorted(old|INPUTS|extra),indent=2)+'\n')
def model_base(stage,geom):
 d=read(f'results/cardiac_subtype_niches/{stage}_{geom}_CM_subtype_neighborhoods.tsv.gz',index_col='cell_id')
 c=read(f'results/cardiac_spatial_reconstruction/{stage}/CM_outcomes_and_anatomy.tsv.gz',index_col=0)
 for k in ['metabolic','sarcomere','excitation_contraction','cycling']:
  d['native_Z__'+k]=c.loc[d.index,'native_Z__'+k]
 d['core_Z__cycling']=-d.core_Z__withdrawal
 extra=read(OUT/(stage+'_CM_Myl2_sensitivity.tsv.gz'),index_col='cell_id')
 return d.join(extra)
def preparation():
 ann,x,raw,genes=refs();mean=(x/np.maximum(x.sum(1),1)[:,None]).mean(0)
 assert set(sum(PANELS.values(),[])).issubset(genes)
 for p in ['results/fibroblast_Xenium_reviewed_refit/RNA_clean_fibroblast_annotations.tsv.gz','results/xenium_ventricular_mapping/reference_shared_panel_counts.npz']:INPUTS.add(str(ROOT/p))
 coverage=[];source_summaries=[];source_records=[];summary=[]
 cal=read('results/cardiac_spatial_reconstruction/score_reference_calibration.tsv');cal=cal[cal.view.eq('core')&cal.scope.eq('all')].set_index('gene')
 for stage in STAGES:
  print('REFIT_PREPARE',stage,flush=True)
  fb,q=queries(stage,genes)
  labels=read(f'results/fibroblast_multimethod_identity/{stage}_multimethod_FB_annotations.tsv.gz',index_col='cell_id')
  f=read(f'results/cardiac_spatial_pipeline/{stage}/atlas_annotations.tsv.gz',index_col=0)
  assert labels.index.equals(fb.index)
  eligible=labels.program_source_eligible.to_numpy(bool)
  present=set(ann.loc[ann.age.eq(stage),'fibroblast_family'].str.split('-').str[0])
  for variant,col in LABELS.items():
   for fam in FAMILIES:
    assessed=fam in present;n=int((eligible&labels[col].eq(fam+'-associated').to_numpy()).sum())
    coverage.append({'stage':stage,'label_layer':variant,'family':fam,'eligible_FB_profiles':int(eligible.sum()),'reference_assessable':assessed,
      'candidate_profiles':n if assessed else np.nan,'table_display':str(n) if assessed else 'Not assessable',
      'meaning':'Assignment coverage, not tissue abundance or independently measured accuracy'})
  # Matched-depth random loss tests information loss without changing boundaries.
  rng=np.random.default_rng(SEED+int(stage[1:]));ratio=np.minimum(1,q['core'].sum(1)/np.maximum(q['refined'].sum(1),1))
  thinned=[rng.binomial(q['refined'].astype(np.int64),ratio[:,None]).astype(float) for _ in range(3)]
  counts={k:q[k] for k in ['refined','core','native']};counts['depthmatched']=thinned[0]
  scores={};info={};depth={};detected={}
  for view,xx in counts.items():
   m=score_modules(xx,genes,PANELS,mean)
   if view=='depthmatched':
    reps=[m]+[score_modules(t,genes,PANELS,mean) for t in thinned[1:]]
    m=sum(reps)/3
   else:
    old=read(f'results/fibroblast_multimethod_identity/{stage}_{view}_published_module_scores.tsv.gz',index_col='cell_id')
    for mod in MODULES[:4]:assert np.allclose(m[mod+'__rank_excess'],old.loc[labels.index,mod+'__rank_excess'],atol=1e-7)
   depth[view]=xx.sum(1);detected[view]=(xx>0).sum(1)
   if view=='depthmatched':
    depth[view]=np.mean([t.sum(1) for t in thinned],axis=0);detected[view]=np.mean([(t>0).sum(1) for t in thinned],axis=0)
   info[view]=eligible&(depth[view]>=20)&(detected[view]>=8)
   # Core and depth-thinned views use identical informative-source masks.
   if view=='depthmatched':info[view]=info['core']
   scores[view]=np.column_stack([m[k+'__rank_excess'].to_numpy() for k in MODULES])
   for j,mod in enumerate(MODULES):
    markers=m[mod+'__measured_markers'].to_numpy();flag=(m[mod+'__control_percentile'].to_numpy()>=.9)&(markers>=2)
    source_summaries.append({'stage':stage,'view':view,'program':mod,'eligible_FB_profiles':int(eligible.sum()),'informative_FB_profiles':int(info[view].sum()),
       'informative_with_two_markers':int((info[view]&(markers>=2)).sum()),'enriched_profiles':int((info[view]&flag).sum()),
       'mean_rank_excess':float(scores[view][info[view],j].mean()),'genes':';'.join(PANELS[mod])})
   record=pd.DataFrame(scores[view],index=labels.index,columns=['program__'+k for k in MODULES]);record['informative']=info[view]
   record['view_depth']=depth[view];record['view_detected_genes']=detected[view];record['cell_id']=labels.index
   record['stage']=stage;record['view']=view;source_records.append(record)
  # Remove Myl2 from both total and structural outcomes without changing selection.
  path=ROOT/f'results/xenium_nuclear_anchored_segmentation/{stage}/nuclear_anchored_counts.h5ad';INPUTS.add(str(path))
  b=ad.read_h5ad(path,backed='r');cm=read(f'results/cardiac_subtype_niches/{stage}_proximity50_CM_subtype_neighborhoods.tsv.gz',index_col='cell_id')
  ix=b.obs_names.get_indexer(cm.index);gi=b.var_names.get_loc('Myl2');yy=b.layers['high_confidence'][ix][:,gi].toarray().ravel();b.file.close()
  zz=(np.log1p(1e4*yy/np.maximum(cm.core_panel_depth.to_numpy(),1))-cal.loc['Myl2','reference_mean'])/cal.loc['Myl2','reference_SD']
  save(pd.DataFrame({'cell_id':cm.index,'core_Myl2_gene_Z':zz,'core_Z__maturation20_noMyl2':(21*cm.core_Z__maturation21.to_numpy()-zz)/20,
        'core_Z__sarcomere_noMyl2':(3*cm.core_Z__sarcomere.to_numpy()-zz)/2}),stage+'_CM_Myl2_sensitivity.tsv.gz')
  for geom in GEOMETRIES:
   base=read(f'results/cardiac_subtype_niches/{stage}_{geom}_CM_subtype_neighborhoods.tsv.gz',index_col='cell_id')
   path=ROOT/f'results/cardiac_spatial_pipeline/{stage}/{geom}_graph.npz';INPUTS.add(str(path));g=sparse.load_npz(path).tocsr()[f.index.get_indexer(base.index)]
   gb=g[:,f.index.get_indexer(labels.index)]
   arrays={};arrays['eligible_FB_neighbors']=np.asarray(gb@eligible.astype(float)).ravel()
   assert np.allclose(arrays['eligible_FB_neighbors'],base.eligible_FB_neighbors)
   for variant,col in LABELS.items():
    for fam in FAMILIES:
     mask=eligible&labels[col].eq(fam+'-associated').to_numpy()
     arrays[variant+'__'+fam]=np.asarray(gb@mask.astype(float)).ravel()/np.maximum(arrays['eligible_FB_neighbors'],1)
   for view in counts:
    good=info[view];den=np.asarray(gb@good.astype(float)).ravel()
    arrays[view+'__information_fraction']=den/np.maximum(arrays['eligible_FB_neighbors'],1)
    arrays[view+'__informative_neighbors']=den
    for k,val in [('log_depth',np.log1p(depth[view])),('detected_genes',detected[view])]:arrays[view+'__'+k]=np.asarray(gb@(good*val)).ravel()/np.maximum(den,1)
    sums=gb@(scores[view]*good[:,None])
    for j,mod in enumerate(MODULES):arrays[view+'__'+mod]=sums[:,j]/np.maximum(den,1)
   np.savez_compressed(OUT/(stage+'_'+geom+'_neighborhood_exposures.npz'),cell_ids=base.index.to_numpy(str),**arrays)
   summary.append({'stage':stage,'geometry':geom,'CM_profiles':len(base),'base_eligible_CM':int(base.analysis_eligible.sum()),
      'common_informative_FB_CM':int((base.analysis_eligible.to_numpy()&(arrays['refined__informative_neighbors']>=1)&(arrays['core__informative_neighbors']>=1)&(arrays['native__informative_neighbors']>=1)).sum())})
 save(pd.DataFrame(coverage),'fibroblast_assignment_coverage.tsv');save(pd.DataFrame(source_summaries),'fibroblast_program_coverage.tsv')
 save(pd.concat(source_records,ignore_index=True),'fibroblast_per_profile_programs.tsv.gz');save(pd.DataFrame(summary),'CM_neighborhood_cohort_counts.tsv')
 save(pd.DataFrame([{'program':k,'genes':';'.join(v),'measurement':'Measured rank excess over 200 reference-abundance-matched control sets','functional_measurement':False} for k,v in PANELS.items()]),'fibroblast_program_definitions.tsv')
 export_inputs()

def nuisance(d,adjustment,view=None):
 cov={'log_CM_capture':np.log1p(d.CM_outcome_excluded_depth),'log_CM_area':np.log1p(d.cell_area),'nucleus_missing':d.nucleus_area.isna().astype(float),
      'log_neighbors':np.log1p(d.neighbors_resolved),'log_FB_count':np.log1p(d.eligible_FB_neighbors),'log_EC_count':np.log1p(d.eligible_EC_neighbors),
      'FB_capture':d.FB_mean_log_capture,'EC_capture':d.EC_mean_log_capture,**{k:d[k] for k in d if k.startswith('broad__')}}
 if view:
  for key in ['information_fraction','log_depth','detected_genes']:cov['FB_'+key]=d[view+'__'+key]
 if adjustment!='technical':
  for k in ['arterial_distance_um','endocardial_proxy_distance_um','distance_to_large_void_um','distance_to_outer_um']:cov['log_'+k]=np.log1p(d[k])
  for key in ['EC__Capillary-like','EC__Endocardial-like']:cov[key]=d[key]
  for key in sorted(d.morphology_region.dropna().unique())[1:]:cov['morphology__'+key]=d.morphology_region.eq(key).astype(float)
 return pd.DataFrame(cov,index=d.index)

def fit_target(d,target,kind,layer,adjustment,geom,extra_cov=None):
 view=layer if kind=='program' else None
 n=nuisance(d,adjustment,view)
 if kind=='vascular':
  n=n.drop(columns=[target],errors='ignore')
 for col in extra_cov or []:n[col]=d[col]
 y=d[PRIMARY_Y+SENS_Y].to_numpy(float);e=d[target].to_numpy(float)
 # One stage-level scale per run preserves CM-subtype comparisons.
 yy=(y-np.nanmean(y,axis=0))/np.maximum(np.nanstd(y,axis=0),1e-6)
 ee=(e-np.nanmean(e))/max(np.nanstd(e),1e-6)
 nx=n.to_numpy(float);nx=(nx-np.nanmean(nx,axis=0))/np.maximum(np.nanstd(nx,axis=0),1e-6)
 finite=np.isfinite(nx).all(1)&np.isfinite(yy).all(1)&np.isfinite(ee)
 allblocks=sorted(d.loc[finite,'block400'].unique());nb=len(allblocks);blockmap={b:i for i,b in enumerate(allblocks)}
 if nb<8:return [],[],[]
 rng=np.random.default_rng(SEED+int(d.stage.iloc[0][1:]));weights=np.stack([np.bincount(rng.integers(nb,size=nb),minlength=nb) for _ in range(DRAWS)])
 rows=[];fits={};diagnostics=[]
 for state in prior.CM:
  take=finite&d.CM_state.eq(state).to_numpy();z=d.loc[take]
  if len(z)<200 or z.block400.nunique()<8:continue
  if kind!='program' and (e[take]>0).sum()<50:continue
  groups=z.block200.to_numpy() if adjustment=='within_region' else np.full(len(z),'intercept')
  cov=prior.center(nx[take],groups);ex=prior.center(ee[take],groups).ravel();response=prior.center(yy[take],groups)
  kept=np.flatnonzero(cov.std(0)>.015)
  if len(kept):
   _,r,piv=linalg.qr(cov[:,kept],mode='economic',pivoting=True);rank=int((abs(np.diag(r))>abs(r[0,0])*1e-8).sum());kept=kept[piv[:rank]]
  cov=cov[:,kept];resid=ex-cov@np.linalg.lstsq(cov,ex,rcond=1e-8)[0]
  if resid.std()<.02:continue
  design=np.column_stack([cov,ex]);p=design.shape[1];q=response.shape[1]
  xx=np.zeros((nb,p,p));xy=np.zeros((nb,p,q))
  codes=np.array([blockmap[k] for k in z.block400]);
  for code in np.unique(codes):
   ix=codes==code;xx[code]=design[ix].T@design[ix];xy[code]=design[ix].T@response[ix]
  point=np.linalg.pinv(xx.sum(0),rcond=1e-9)@xy.sum(0)
  A=np.einsum('bn,nij->bij',weights,xx);B=np.einsum('bn,nij->bij',weights,xy)
  boot=np.linalg.pinv(A,rcond=1e-8)@B
  beta=point[-1];boots=boot[:,-1];fits[state]=(beta,boots)
  diagnostics.append({'stage':d.stage.iloc[0],'geometry':geom,'source_kind':kind,'source_layer':layer,'adjustment':adjustment,'target':target,'CM_state':state,
    'CM_profiles':len(z),'blocks400':z.block400.nunique(),'design_rank':p,'design_columns_before_QR':len(n.columns)+1,'condition_number':float(np.linalg.cond(design)),
    'target_residual_SD':resid.std(),'local_regions':z.block200.nunique(),'nonzero_exposure_CM':int((e[take]!=0).sum())})
  for j,outcome in enumerate(PRIMARY_Y+SENS_Y):
   low,high=np.quantile(boots[:,j],[.025,.975]);se=boots[:,j].std(ddof=1)
   rows.append({'stage':d.stage.iloc[0],'geometry':geom,'source_kind':kind,'source_layer':layer,'adjustment':adjustment,'target':target,'CM_state':state,'outcome':outcome,
      'coefficient':beta[j],'CI_low':low,'CI_high':high,'within_section_p':2*norm.sf(abs(beta[j])/max(se,1e-9)),
      'CM_profiles':len(z),'blocks400':z.block400.nunique(),'primary_outcome':outcome in PRIMARY_Y,'scope':'Within-section observational association; one section per age; standardized within-stage run units'})
 dif=[]
 if all(k in fits for k in ['v.CM','Ankrd1.CM']):
  a,ab=fits['Ankrd1.CM'];v,vb=fits['v.CM'];delta=a-v;boot=ab-vb
  for j,outcome in enumerate(PRIMARY_Y+SENS_Y):
   lo,hi=np.quantile(boot[:,j],[.025,.975]);dif.append({'stage':d.stage.iloc[0],'geometry':geom,'source_kind':kind,'source_layer':layer,'adjustment':adjustment,'target':target,'outcome':outcome,
    'Ankrd1_minus_working':delta[j],'CI_low':lo,'CI_high':hi,'scope':'Shared spatial-block resamples; not independent-heart interaction validation'})
 return rows,dif,diagnostics

def bh(p):
 p=np.asarray(p);order=np.argsort(p);out=np.empty(len(p));out[order]=np.minimum(1,np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]);return out

def fitting():
 coverage=read(OUT/'fibroblast_assignment_coverage.tsv');pcov=read(OUT/'fibroblast_program_coverage.tsv')
 rows=[];diff=[];diag=[]
 for stage in STAGES:
  for geom in GEOMETRIES:
   d=model_base(stage,geom)
   with np.load(OUT/(stage+'_'+geom+'_neighborhood_exposures.npz')) as b:
    assert np.array_equal(b['cell_ids'],d.index.to_numpy(str))
    exposures=pd.DataFrame({k:b[k] for k in b.files if k!='cell_ids'},index=d.index)
    d=pd.concat([d.drop(columns=d.columns.intersection(exposures.columns)),exposures],axis=1)
   # Same CMs in native/refined/core/thinned comparisons within each geometry.
   common=d.analysis_eligible&(d.refined__informative_neighbors>=1)&(d.core__informative_neighbors>=1)&(d.native__informative_neighbors>=1)
   d=d.loc[common].copy()
   variants=list(LABELS) if geom=='proximity50' else ['candidate','crosschecked']
   for layer in variants:
    cov=coverage[coverage.stage.eq(stage)&coverage.label_layer.eq(layer)&coverage.reference_assessable&coverage.candidate_profiles.ge(30)]
    for fam in cov.family:
     target=layer+'__'+fam
     for adjust in (['technical','anatomy','within_region'] if geom=='proximity50' and layer=='candidate' else ['within_region']):
      a,b,c=fit_target(d,target,'identity',layer,adjust,geom);rows+=a;diff+=b;diag+=c
   for view in (['refined','core','native','depthmatched'] if geom=='proximity50' else ['refined','core']):
    mods=pcov.loc[pcov.stage.eq(stage)&pcov.view.eq(view)&pcov.informative_with_two_markers.ge(30),'program']
    for mod in mods:
     target=view+'__'+mod
     for adjust in (['technical','anatomy','within_region'] if geom=='proximity50' and view=='refined' else ['within_region']):
      a,b,c=fit_target(d,target,'program',view,adjust,geom);rows+=a;diff+=b;diag+=c
   # Reassess the endocardial/capillary relationships after updated FB adjustment.
   for layer in (['candidate','programs'] if geom=='proximity50' else ['candidate']):
    fb_cov=[layer+'__'+f for f in FAMILIES if d.get(layer+'__'+f,pd.Series([0])).std()>1e-6] if layer=='candidate' else ['refined__'+m for m in MODULES]
    for target in ['EC__Capillary-like','EC__Endocardial-like']:
     a,b,c=fit_target(d,target,'vascular',layer,'within_region',geom,extra_cov=fb_cov);rows+=a;diff+=b;diag+=c
   print('REFIT_MODELS',stage,geom,'cumulative_coefficients',len(rows),flush=True)
 r=pd.DataFrame(rows)
 r['primary_analysis']=r.geometry.eq('proximity50')&r.adjustment.eq('within_region')&r.primary_outcome&(
  (r.source_kind.eq('identity')&r.source_layer.eq('candidate'))|(r.source_kind.eq('program')&r.source_layer.eq('refined'))|(r.source_kind.eq('vascular')&r.source_layer.eq('candidate')))
 r['primary_grid_BH_q']=np.nan;r.loc[r.primary_analysis,'primary_grid_BH_q']=bh(r.loc[r.primary_analysis,'within_section_p'])
 save(r,'CM_fibroblast_refitted_coefficients.tsv');save(pd.DataFrame(diff),'Ankrd1_working_slope_comparisons.tsv');save(pd.DataFrame(diag),'model_design_diagnostics.tsv')
 export_inputs()

def transport_descriptors():
 # Keep frozen molecular distances/plans; update held-out FB descriptions only.
 old=read('results/cardiac_subtype_niches/CM_identity_restricted_OT_neighborhoods.tsv')
 clouds=read('results/cardiac_niche_dynamics/Wormhole_heldout_clouds.tsv');sources=read(OUT/'fibroblast_per_profile_programs.tsv.gz')
 rows=[]
 for stage in STAGES:
  lab=read(f'results/fibroblast_multimethod_identity/{stage}_multimethod_FB_annotations.tsv.gz',index_col='cell_id')
  with np.load(ROOT/f'results/cardiac_niche_dynamics/{stage}/nonCM_state_points.npz') as b:ids=b['source_ids'].astype(str)
  path=ROOT/f'results/cardiac_niche_dynamics/{stage}/state_cloud_membership.npz';INPUTS.add(str(path));member=sparse.load_npz(path).tocsr()
  prog=sources[sources.stage.eq(stage)&sources.view.eq('refined')].set_index('cell_id')
  for _,r in old[old.stage.eq(stage)].iterrows():
   focal=int(clouds.loc[int(r.cloud_index),'focal_index']);ix=member.getrow(focal).indices
   if len(ix)>32:ix=np.random.default_rng(947+focal).choice(ix,32,replace=False)
   fbids=pd.Index(ids[ix]).intersection(lab.index);eligibleids=fbids[lab.loc[fbids,'program_source_eligible'].to_numpy(bool)]
   record={'stage':stage,'cell_id':r.cell_id,'cloud_index':int(r.cloud_index),'CM_state':r.CM_state,'distance_stratum':r.distance_stratum,
    'sampled_sources':len(ix),'eligible_FB_sources':len(eligibleids),'same_frozen_neighborhood':True}
   for fam in FAMILIES:record['candidate__'+fam]=float(lab.loc[eligibleids,'recommended_candidate'].eq(fam+'-associated').mean()) if len(eligibleids) else np.nan
   good=eligibleids[prog.loc[eligibleids,'informative'].to_numpy(bool)]
   record['informative_FB_sources']=len(good)
   for mod in MODULES:record['program__'+mod]=prog.loc[good,'program__'+mod].mean() if len(good) else np.nan
   rows.append(record)
 out=pd.DataFrame(rows);assert len(out)==96;save(out,'updated_frozen_transport_neighborhood_descriptors.tsv')
 # Descriptor differences can be compared after matching, without training on them.
 plans=read('results/cardiac_subtype_niches/CM_identity_restricted_OT_correspondence.tsv');o=out.set_index('cloud_index');summary=[]
 for key,g in plans.groupby(['stage_a','stage_b','CM_state','distance_stratum']):
  for mod in MODULES:
   a=o.loc[g.cloud_a,'program__'+mod].to_numpy();b=o.loc[g.cloud_b,'program__'+mod].to_numpy();good=np.isfinite(a)&np.isfinite(b)
   for method,col in [('exact','exact_correspondence_mass'),('Wormhole','Wormhole_correspondence_mass')]:
    mass=g[col].to_numpy();supported=float(mass[good].sum())
    summary.append({'stage_a':key[0],'stage_b':key[1],'CM_state':key[2],'distance_stratum':key[3],'program':mod,'method':method,
     'supported_transport_mass':supported,'mean_target_minus_source':float(np.sum(mass[good]*(b[good]-a[good]))/supported) if supported else np.nan,
     'scope':'Updated held-out descriptor on frozen matches; not ancestry, functional change or new encoder training'})
 save(pd.DataFrame(summary),'updated_transport_program_comparisons.tsv');export_inputs()

def main():
 p=argparse.ArgumentParser();p.add_argument('--step',choices=['prepare','models','transport','all'],default='all');a=p.parse_args();OUT.mkdir(parents=True,exist_ok=True);t=time.time()
 if a.step in ['prepare','all']:preparation()
 if a.step in ['models','all']:fitting()
 if a.step in ['transport','all']:transport_descriptors()
 (OUT/('execution_'+a.step+'.json')).write_text(json.dumps({'completed':True,'seconds':time.time()-t,'seed':SEED,'bootstrap_draws':DRAWS,'new_biological_replicates':0,'new_segmentation':False},indent=2)+'\n')
 print('REFIT_COMPLETE',a.step,time.time()-t,flush=True)
if __name__=='__main__':main()
