"""Retest vascular proximity before conditioning on possible molecular mediators.

Frozen cell identities/segmentation; no oxygen imputation or causal mediation.
All intervals describe one observed section per age, not animal replication.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','2')
os.environ.setdefault('MPLCONFIGDIR','/tmp/capillary-retest-mpl')
from pathlib import Path
import json, hashlib
import numpy as np
import pandas as pd
import h5py
from anndata.io import read_elem,sparse_dataset
from scipy import sparse,linalg
from scipy.spatial import cKDTree
from scipy.stats import t,spearmanr,pearsonr
from sklearn.preprocessing import SplineTransformer
import cardiac_subtype_niches as prior
from assess_CM_maturation21_neighborhoods import PANEL,COMPONENTS,bh

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/capillary_maturation_retest';OUT.mkdir(parents=True,exist_ok=True)
INPUTS=set();STAGES=['P0','P7','P14','P21']
def read(p,**kw):
 p=ROOT/p;INPUTS.add(str(p));return pd.read_csv(p,sep='\t',**kw)
def save(d,p,index=False):d.to_csv(OUT/p,sep='\t',index=index,float_format='%.8g',na_rep='NA')
def std(x):return (x-np.mean(x,axis=0))/np.maximum(np.std(x,axis=0),1e-9)
def export_cm(cm,stage):
 keep=['stage_cell_id','x_centroid','y_centroid','author_atlas_cell_state','CM_state','strict','capture',
  'cell_area','nucleus_missing','core_panel_depth','native_panel_depth','core_depth50','anatomy_block_200um',
  'morphology_region','distance_to_large_void_um','distance_to_outer_um']
 keep += [v+'_'+s for v in ['core','native'] for s in ['maturation21','metabolic','sarcomere','excitation_contraction','cycling','matched3_100__maturation21']]
 keep += [k for k in cm if k.startswith(('near_','distance_public','distance_QC','apposition_','core_contact_','native_contact_'))]
 save(cm[list(dict.fromkeys(keep))],stage+'_CM_retest_values.tsv.gz',True)

def fit(frame,exposures,outcomes,adjustment,block=400):
 """OLS with pivoted rank check and finite-sample clustered t uncertainty."""
 x=frame[exposures].to_numpy(float);y=frame[outcomes].to_numpy(float)
 cov=[];base=adjustment.replace('_spline','')
 if base!='unadjusted':
  cov=[np.log1p(frame.capture),np.log1p(frame.cell_area),frame.nucleus_missing.astype(float)]
  if adjustment.endswith('_spline'):
   cov=[SplineTransformer(n_knots=4,degree=3,include_bias=False).fit_transform(np.column_stack(cov[:2])),cov[2]]
 if base in ['state_anatomy','within_region_state']:
  cov += [np.log1p(frame.distance_to_large_void_um),np.log1p(frame.distance_to_outer_um)]
  cov += [pd.get_dummies(frame.CM_state,dtype=float).to_numpy(),pd.get_dummies(frame.morphology_region,dtype=float).to_numpy()]
 c=np.column_stack(cov) if cov else np.zeros((len(frame),0))
 valid=np.isfinite(x).all(1)&np.isfinite(y).all(1)&np.isfinite(c).all(1)
 z=frame.loc[valid];x=std(x[valid]);y=std(y[valid]);c=std(c[valid]) if c.shape[1] else c[valid]
 if base=='within_region_state':
  group=z.anatomy_block_200um.astype(str)+'::'+z.CM_state.astype(str)
  x=prior.center(x,group);y=prior.center(y,group);c=prior.center(c,group)
 else:c=np.column_stack([np.ones(len(z)),c])
 if c.shape[1]:
  _,r,piv=linalg.qr(c,mode='economic',pivoting=True)
  rank=np.sum(abs(np.diag(r))>max(abs(r[0,0]),1e-12)*1e-9)
  c=c[:,piv[:rank]]
 design=np.column_stack([c,x]);assert np.linalg.matrix_rank(design)==design.shape[1]
 bread=np.linalg.pinv(design.T@design,rcond=1e-10);beta=bread@design.T@y;resid=y-design@beta
 groups=(z.x_centroid//block).astype(int).astype(str)+'_'+(z.y_centroid//block).astype(int).astype(str)
 codes=pd.factorize(groups)[0];ng=codes.max()+1
 assert ng>=8
 # Cluster score contributions; correlations between outcomes do not affect per-outcome intervals.
 scores=np.zeros((ng,design.shape[1],len(outcomes)))
 for k in range(ng):
  ix=codes==k;scores[k]=design[ix].T@resid[ix]
 correction=ng/(ng-1)*(len(z)-1)/(len(z)-design.shape[1])
 rows=[]
 for i,e in enumerate(exposures):
  j=c.shape[1]+i;v=np.einsum('p,gpq->gq',bread[j],scores)
  se=np.sqrt(np.sum(v*v,axis=0)*correction);b=beta[j];ci=t.ppf(.975,ng-1)*se
  for k,o in enumerate(outcomes):
   rows.append(dict(exposure=e,outcome=o,adjustment=adjustment,block_um=block,CM_profiles=len(z),spatial_blocks=ng,
    beta_SD=b[k],CI_low=b[k]-ci[k],CI_high=b[k]+ci[k],within_section_p=2*t.sf(abs(b[k]/max(se[k],1e-15)),ng-1),
    design_rank=design.shape[1],condition_number=np.linalg.cond(design)))
 return rows

def spline_sensitivity():
 rows=[]
 for stage in STAGES:
  cm=pd.read_csv(OUT/(stage+'_CM_retest_values.tsv.gz'),sep='\t',index_col=0)
  sub=cm[cm.strict]
  for exposure in ['near_public_capillary','near_QC_capillary','apposition_author_capillary','apposition_QC_capillary']:
   for adjustment in ['state_anatomy_spline','within_region_state_spline']:
    for block in [400,800]:
     rows += [dict(stage=stage,cohort='strict_ventricular',**r) for r in fit(sub,[exposure],['core_maturation21','native_maturation21'],adjustment,block)]
  deep=sub[sub.core_depth50.notna()]
  rows += [dict(stage=stage,cohort='strict_depth50',**r) for r in fit(deep,['near_QC_capillary'],['core_maturation21','core_depth50'],'within_region_state_spline')]
 save(pd.DataFrame(rows),'flexible_capture_sensitivity.tsv')

def contact_spline_sensitivity():
 rows=[];checks=[]
 for stage in STAGES:
  f,d=prior.stage_data(stage);cm=pd.read_csv(OUT/(stage+'_CM_retest_values.tsv.gz'),sep='\t',index_col=0)
  cap=(f.atlas_identity.eq('Endothelial')&f.program_source_eligible&f.EC_nuclear_type.eq('Capillary-like')).to_numpy(float)
  cms=(f.atlas_identity.eq('Cardiomyocyte')&f.pass_expression_QC).to_numpy(float)
  si=f.index.get_indexer(d.index);G=sparse.load_npz(ROOT/f'results/cardiac_spatial_pipeline/{stage}/contact_graph.npz').tocsr()[si]
  for view in ['core','native']:
   path=ROOT/(f'results/xenium_nuclear_anchored_segmentation/{stage}/nuclear_anchored_counts.h5ad' if view=='core' else f'results/xenium_ventricular_mapping/{stage}/{stage}_xenium_counts_and_candidate_mapping.h5ad')
   with h5py.File(path) as h:
    ids=read_elem(h['obs']).index;genes=read_elem(h['var']).index
    x=sparse_dataset(h['layers/high_confidence'] if view=='core' else h['X'])[ids.get_indexer(f.index)].tocsr()
    depth=np.asarray(x[:,~genes.isin(set(PANEL)|set(COMPONENTS['cycling']))].sum(1)).ravel()
    logpt=np.log1p(x[:,genes.get_loc('Ptprm')].toarray().ravel()*1e4/np.maximum(depth,1));del x
   m=d.copy();m['CM_state']=m.author_atlas_cell_state;m['capture']=depth[si]
   m['anatomy_block_200um']=cm.loc[m.index,'anatomy_block_200um'];m['nucleus_missing']=m.nucleus_area.isna()
   n=np.asarray(G@cap).ravel();nc=np.asarray(G@cms).ravel()
   m['near_capillary']=cm.loc[m.index,'near_QC_capillary'];m['cap_count']=np.log1p(n);m['CM_neighbor_count']=np.log1p(nc)
   m['cap_capture']=np.asarray(G@(cap*np.log1p(depth))).ravel()/np.maximum(n,1)
   m['cap_Ptprm']=np.asarray(G@(cap*logpt)).ravel()/np.maximum(n,1)
   m['receiver_Ptprm']=logpt[si]
   m['CM_neighbor_Ptprm']=np.asarray(G@(cms*logpt)).ravel()/np.maximum(nc,1)
   m['cap_CM_Ptprm_product']=m.cap_Ptprm*m.receiver_Ptprm
   exposures=['near_capillary','cap_count','cap_capture','CM_neighbor_count','CM_neighbor_Ptprm','receiver_Ptprm','cap_Ptprm','cap_CM_Ptprm_product']
   yy=['Y_'+s for s in ['maturation21','metabolic','sarcomere','excitation_contraction','cycling']]
   for key in yy:m[key]=m[view+'_Z__'+key[2:]]
   for block in [400,800]:
    rows += [dict(stage=stage,view=view,geometry='contact',**r) for r in fit(m,exposures,yy,'within_region_state_spline',block)]
   # Rebuild the earlier model from original count files as a consistency check.
   old=pd.read_csv(OUT/'PTPRM_and_Vegfa_regressions.tsv',sep='\t')
   old=old[old.stage.eq(stage)&old.view.eq(view)&old.geometry.eq('contact')&old.block_um.eq(400)].set_index(['exposure','outcome'])
   new=pd.DataFrame(fit(m,exposures,yy,'within_region_state')).set_index(['exposure','outcome'])
   error=float(abs(old.beta_SD-new.beta_SD).max());assert error<1e-7
   checks.append(dict(stage=stage,view=view,original_count_rebuilt_model_max_abs_error=error,status='PASS'))
 save(pd.DataFrame(rows),'PTPRM_flexible_capture_sensitivity.tsv')
 save(pd.DataFrame(checks),'PTPRM_count_reconstruction_checks.tsv')

def main():
 prox=[];mechanism=[];coverage=[];curves=[];graphrows=[];matched=[];rna=[]
 cal=read('results/cardiac_spatial_reconstruction/score_reference_calibration.tsv').set_index(['view','scope','gene']).sort_index()
 for stage in STAGES:
  print('RETEST',stage,flush=True)
  INPUTS.add(str(ROOT/f'results/cardiac_spatial_pipeline/{stage}/atlas_annotations.tsv.gz'))
  f,d=prior.stage_data(stage)
  c=read(f'results/cardiac_spatial_reconstruction/{stage}/CM_outcomes_and_anatomy.tsv.gz',index_col=0)
  q=read(f'results/cardiac_spatial_reconstruction/{stage}/all_profile_observations.tsv.gz',index_col=0)
  rings=read(f'results/cardiac_spatial_reconstruction/{stage}/reconstructed_neighbor_levels.tsv.gz',index_col=0)
  scores=read(f'results/CM_maturation21_nonCM_context/{stage}/CM_maturation21_scores.tsv.gz',index_col=0)
  boundary=read(f'results/cardiac_spatial_reconstruction/{stage}/vascular_interface_and_matrix_context.tsv.gz',index_col=0)
  cm=c[c.reconstructed_Figure4e_CM&c.core_panel_depth.ge(20)].copy()
  cm['strict']=cm.index.isin(d.index);cm['CM_state']=cm.author_atlas_cell_state
  cm['capture']=f.loc[cm.index,'CM_outcome_excluded_depth'];cm['nucleus_missing']=cm.nucleus_area_missing
  cm['core_depth50']=scores.loc[cm.index,'core_Z__depth50_mean3']
  for view in ['core','native']:
   for score in ['metabolic','sarcomere','excitation_contraction','cycling','maturation21']:cm[view+'_'+score]=cm[view+'_Z__'+score]
  masks={
   'public_capillary':f.author_atlas_cell_state.eq('cap.EC'),
   'QC_capillary':f.atlas_identity.eq('Endothelial')&f.program_source_eligible&f.EC_nuclear_type.eq('Capillary-like')}
  xy=f[['x_centroid','y_centroid']].to_numpy()
  for name,mask in masks.items():
   dist=cKDTree(xy[mask]).query(cm[['x_centroid','y_centroid']])[0]
   cm['distance_'+name]=dist;cm['near_'+name]=-np.log1p(dist)
   coverage.append(dict(stage=stage,source_definition=name,source_profiles=int(mask.sum()),broad_CMs=len(cm),strict_CMs=int(cm.strict.sum())))
  for name in ['author_capillary','QC_capillary']:
   cm['apposition_'+name]=boundary.loc[cm.index,name+'__apposition_fraction_2um']
  for cohort,sub in [('reconstructed_paper_labels',cm),('strict_ventricular',cm[cm.strict])]:
   outcomes=['core_maturation21','native_maturation21','core_metabolic','core_sarcomere','core_excitation_contraction','core_cycling',
    'core_matched3_100__maturation21','native_matched3_100__maturation21']
   for exposure in ['near_public_capillary','near_QC_capillary']:
    for adjust in ['unadjusted','technical','state_anatomy','within_region_state']:
     for block in [400,800]:
      prox += [dict(stage=stage,cohort=cohort,**r) for r in fit(sub,[exposure],outcomes,adjust,block)]
   for exposure in ['apposition_author_capillary','apposition_QC_capillary']:
    prox += [dict(stage=stage,cohort=cohort,**r) for r in fit(sub,[exposure],outcomes,'within_region_state')]
   # Descriptive bins remain at the cell level, with within-bin cell spread.
   for source in masks:
    bins=pd.cut(sub['distance_'+source],[0,10,20,30,50,75,100,200,np.inf],include_lowest=True)
    for label,g in sub.groupby(bins,observed=True):
     curves.append(dict(stage=stage,cohort=cohort,source=source,distance_bin=str(label),CMs=len(g),
      mean_distance=g['distance_'+source].mean(),mean_core21=g.core_maturation21.mean(),mean_native21=g.native_maturation21.mean()))
   for graph in ['paper_like6NN__cap.EC','tissue_Delaunay__cap.EC']:
    level=rings.loc[sub.index,graph];means=sub.core_maturation21.groupby(level).mean()
    rr=pearsonr(means.index,means) if len(means)>2 else (np.nan,np.nan)
    graphrows.append(dict(stage=stage,cohort=cohort,graph=graph,level_mean_Pearson_r=rr[0],
     cell_level_Spearman_r=spearmanr(level,sub.core_maturation21,nan_policy='omit').statistic,
     populated_levels=len(means),CMs=len(sub),interpretation='Descriptive; level means are not independent biological replicates'))
  # Exact same previous matched cells, independently recomputed score contrasts.
  pp=read('results/CM_maturation21_nonCM_context/matched_capillary_pairs_used.tsv.gz')
  pp=pp[pp.stage.eq(stage)]
  for cohort,g in pp.groupby('cohort'):
   a=c.index.get_indexer(g.high_apposition_CM);b=c.index.get_indexer(g.zero_apposition_CM);assert (a>=0).all() and (b>=0).all()
   for view in ['core','native']:
    yy=c[view+'_Z__maturation21'];delta=yy.iloc[a].to_numpy()-yy.iloc[b].to_numpy()
    matched.append(dict(stage=stage,cohort=cohort,view=view,matched_pairs=len(g),raw_score_difference=delta.mean()))
  # Read deposited and boundary-supported counts; same measured genes in both views.
  for view in ['core','native']:
   path=ROOT/(f'results/xenium_nuclear_anchored_segmentation/{stage}/nuclear_anchored_counts.h5ad' if view=='core' else f'results/xenium_ventricular_mapping/{stage}/{stage}_xenium_counts_and_candidate_mapping.h5ad')
   INPUTS.add(str(path))
   with h5py.File(path) as h:
    ids=read_elem(h['obs']).index;genes=read_elem(h['var']).index;ix=ids.get_indexer(f.index);assert (ix>=0).all()
    x=sparse_dataset(h['layers/high_confidence'] if view=='core' else h['X'])[ix].tocsr()
    exclude=set(PANEL)|set(COMPONENTS['cycling']);depth=np.asarray(x[:,~genes.isin(exclude)].sum(1)).ravel()
    pt=x[:,genes.get_loc('Ptprm')].toarray().ravel();ve=x[:,genes.get_loc('Vegfa')].toarray().ravel()
    logpt=np.log1p(pt*1e4/np.maximum(depth,1));logve=np.log1p(ve*1e4/np.maximum(depth,1))
    # Fixed-reference leave-one-gene-out sensitivity of the proximity association itself.
    ci=f.index.get_indexer(cm.index);total=np.asarray(x.sum(1)).ravel();pc=x[ci][:,genes.get_indexer(PANEL)].toarray()
    ref=cal.loc[(view,'all')].reindex(PANEL);z=(np.log1p(pc*1e4/np.maximum(total[ci],1)[:,None])-ref.reference_mean.to_numpy())/ref.reference_SD.to_numpy()
    np.testing.assert_allclose(z.mean(1),cm[view+'_maturation21'],atol=2e-6)
    for j,gene in enumerate(PANEL):cm[view+'_omit_'+gene]=(z.sum(1)-z[:,j])/20
    del x
   si=f.index.get_indexer(d.index);m=d.copy();m['CM_state']=m.author_atlas_cell_state;m['capture']=depth[si]
   m['anatomy_block_200um']=cm.loc[m.index,'anatomy_block_200um']
   m['nucleus_missing']=m.nucleus_area.isna();m['receiver_Ptprm']=logpt[si];m['receiver_Vegfa']=logve[si]
   for outcome in ['maturation21','metabolic','sarcomere','excitation_contraction','cycling']:m['Y_'+outcome]=m[view+'_Z__'+outcome]
   for geometry in ['proximity50','contact']:
    p=ROOT/f'results/cardiac_spatial_pipeline/{stage}/{geometry}_graph.npz';INPUTS.add(str(p));G=sparse.load_npz(p).tocsr()[si]
    cap=masks['QC_capillary'].to_numpy(float)
    n=np.asarray(G@cap).ravel();m['cap_count']=np.log1p(n)
    m['cap_Ptprm']=np.asarray(G@(cap*logpt)).ravel()/np.maximum(n,1)
    m['cap_capture']=np.asarray(G@(cap*np.log1p(depth))).ravel()/np.maximum(n,1)
    cms=(f.atlas_identity.eq('Cardiomyocyte')&f.pass_expression_QC).to_numpy(float)
    nc=np.asarray(G@cms).ravel();m['CM_neighbor_count']=np.log1p(nc)
    m['CM_neighbor_Ptprm']=np.asarray(G@(cms*logpt)).ravel()/np.maximum(nc,1)
    # Explicit product plus both main terms, with source abundance/capture controls.
    m['cap_CM_Ptprm_product']=m.cap_Ptprm*m.receiver_Ptprm
    m['near_capillary']=cm.loc[m.index,'near_QC_capillary']
    exposures=['near_capillary','cap_count','cap_capture','CM_neighbor_count','CM_neighbor_Ptprm','receiver_Ptprm','cap_Ptprm','cap_CM_Ptprm_product']
    for block in [400,800]:
     mechanism += [dict(stage=stage,view=view,geometry=geometry,**r) for r in fit(m,exposures,['Y_'+s for s in ['maturation21','metabolic','sarcomere','excitation_contraction','cycling']],'within_region_state',block)]
    coverage.append(dict(stage=stage,source_definition=geometry+'_'+view,source_profiles=int(cap.sum()),strict_CMs=len(m),
     CMs_with_capillary=int((n>0).sum()),CM_Ptprm_detected=float((pt[si]>0).mean()),
     CMs_with_Ptprm_positive_capillary=int((np.asarray(G@(cap*(pt>0))).ravel()>0).sum()),
     CMs_both_sides_Ptprm=int(((pt[si]>0)&(np.asarray(G@(cap*(pt>0))).ravel()>0)).sum())))
    if geometry=='contact':
     for key in ['cap_count','cap_Ptprm','receiver_Ptprm','cap_CM_Ptprm_product']:cm.loc[m.index,view+'_contact_'+key]=m[key]
   # An association with VEGFA is not an oxygen estimate.
   mechanism += [dict(stage=stage,view=view,geometry='single_Vegfa_not_oxygen',**r) for r in fit(m,['near_capillary'],['receiver_Vegfa'],'within_region_state')]
  # Omission/fixed-depth analyses on identical eligible strict CMs.
  sub=cm[cm.strict].copy();cols=[v+'_omit_'+g for v in ['core','native'] for g in PANEL]
  for source in ['near_public_capillary','near_QC_capillary']:
   prox += [dict(stage=stage,cohort='strict_ventricular',**r) for r in fit(sub,[source],cols,'within_region_state')]
   deep=sub[sub.core_depth50.notna()]
   prox += [dict(stage=stage,cohort='strict_depth50',**r) for r in fit(deep,[source],['core_maturation21','core_depth50'],'within_region_state')]
  export_cm(cm,stage)
  # Existing full-transcriptome, expression-controlled hypoxia scores stay in RNA space.
  p=ROOT/f'results/ventricular_microenvironment/{stage}/cell_annotations_and_scores.csv.gz';INPUTS.add(str(p))
  rr=pd.read_csv(p,index_col=0)
  for sample,g in rr.groupby('sample',observed=True):
   for a,b in [('z_HIF_response','z_Maturation_reference'),('qc_adjusted_z_HIF_response','qc_adjusted_z_Maturation_reference'),('z_HIF_non_glycolytic','z_Maturation_reference')]:
    rna.append(dict(stage=stage,sample=sample,program=a,outcome=b,nuclei=len(g),Spearman_r=spearmanr(g[a],g[b]).statistic,
     interpretation='Within-library RNA-state association, not spatial capillary exposure or oxygen measurement'))
 result=pd.DataFrame(prox)
 primary=(result.cohort.eq('strict_ventricular')&result.outcome.eq('core_maturation21')&result.adjustment.eq('state_anatomy')&result.block_um.eq(400)&result.exposure.str.startswith('near_'))
 result['primary_analysis']=primary;result['primary_BH_q']=np.nan;result.loc[primary,'primary_BH_q']=bh(result.loc[primary,'within_section_p'])
 save(result,'proximity_regressions.tsv')
 mech=pd.DataFrame(mechanism)
 primary=mech.view.eq('core')&mech.geometry.eq('contact')&mech.block_um.eq(400)&mech.exposure.isin(['receiver_Ptprm','cap_Ptprm','cap_CM_Ptprm_product'])
 mech['primary_PTPRM_analysis']=primary;mech['PTPRM_BH_q']=np.nan;mech.loc[primary,'PTPRM_BH_q']=bh(mech.loc[primary,'within_section_p'])
 save(mech,'PTPRM_and_Vegfa_regressions.tsv');save(pd.DataFrame(coverage),'source_and_detection_coverage.tsv')
 save(pd.DataFrame(curves),'distance_bin_descriptions.tsv');save(pd.DataFrame(graphrows),'graph_ring_descriptions.tsv')
 save(pd.DataFrame(matched),'independently_recomputed_matched_differences.tsv');save(pd.DataFrame(rna),'RNA_hypoxia_state_correlations.tsv')
 support=read('results/joint_RNA_Xenium_signaling/RNA_pair_source_support.tsv.gz')
 save(support[support.pair_id.eq('Ptprm--Ptprm')],'RNA_PTPRM_source_support.tsv')
 programs=json.loads((ROOT/'results/ventricular_microenvironment/program_definitions.json').read_text())['programs']
 measured=set(cal.index.get_level_values('gene'))
 save(pd.DataFrame([dict(program=k,genes=';'.join(programs[k]['genes']),Xenium_measured=';'.join(g for g in programs[k]['genes'] if g in measured),
  measured_count=len(set(programs[k]['genes'])&measured),total_genes=len(programs[k]['genes'])) for k in ['HIF_response','HIF_non_glycolytic','Oxidative_metabolism']]),'oxygen_program_coverage.tsv')
 (OUT/'methods.json').write_text(json.dumps(dict(primary_proximity_tests=int(result.primary_analysis.sum()),primary_PTPRM_tests=int(primary.sum()),
  score='Frozen pooled gene-Z mean of 21 directly measured genes; exact author AddModuleScore not recovered',
  proximity_primary='8 tests: 4 ages x 2 capillary definitions, strict CMs, technical + identity + anatomy adjustment, 400um clustered uncertainty',
  PTPRM_primary='60 tests: 4 ages x 5 outcomes x 3 Ptprm terms, core RNA, operational contact graph',
  contact='Shared refined raster territory boundary, inherited graph; not verified protein interaction',
  oxygen='No oxygen/perfusion measure or validated spatial hypoxia module; VEGFA alone is insufficient',
  uncertainty='CR1 spatial-cluster sandwich with t(df=blocks-1), 400um primary and 800um sensitivity; one section per age; no animal-level inference',
  mediation_claimed=False,new_segmentation=False),indent=2)+'\n')
 (OUT/'input_paths.json').write_text(json.dumps(sorted(INPUTS),indent=2)+'\n')
 spline_sensitivity()
 contact_spline_sensitivity()
 print(result[result.primary_analysis].to_string(index=False),flush=True)
 print(mech[mech.primary_PTPRM_analysis&mech.PTPRM_BH_q.lt(.05)].to_string(index=False),flush=True)
if __name__=='__main__':main()
