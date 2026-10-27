"""Independent reconstruction of molecular support and spatial prediction checks."""
import json
import h5py
import numpy as np
import pandas as pd
from scipy import sparse
from scipy import ndimage
from scipy.spatial import cKDTree
from sklearn.metrics import r2_score
from anndata.io import sparse_dataset,read_elem
from joint_signaling_common import *

def rna_audit():
 support=read(OUT/'RNA_pair_source_support.tsv.gz');expr=read(OUT/'RNA_complex_expression_by_library.tsv.gz')
 selected=read(OUT/'selected_RNA_supported_pairs.tsv');measured=set(read('results/cardiac_spatial_reconstruction/score_reference_calibration.tsv').gene)
 species=set('_'.join(selected.ligand.tolist()+selected.receptor.tolist()).split('_'))
 assert species<=measured and not species&OUTCOME_GENES
 assert not selected.pair_id.duplicated().any()
 corrected=expr[expr.view.eq('corrected')]
 ok=(corrected.nuclei>=30)&(corrected.mean_log1p_CP10k>=.03)&(corrected.fraction_complex_detected>=.02)
 independently=corrected.assign(pass_=ok).groupby(['stage','source_group','species']).agg(libraries=('sample','nunique'),pass_both=('pass_','all'))
 for row in support.itertuples():
  a=independently.loc[(row.stage,'CM',row.receptor)]
  key=(row.stage,row.source_group,row.ligand)
  b=independently.loc[key] if key in independently.index else None
  sender=b is not None and b.libraries==2 and b.pass_both
  receiver=a.libraries==2 and a.pass_both
  assert bool(row.joint_two_library_support)==bool(sender and receiver),row
 # Reconstruct selected source/library means from the original full-genome
 # corrected matrix, independently of the support-export implementation.
 path='/workspace/heart-atlas/results/cell_identities_subtypes/snrna_qc_annotated_subclustered.h5ad';INPUTS.add(path)
 annotation=read(OUT/'RNA_cell_source_annotations.tsv.gz',index_col=0)
 checked=['Nrg1','Erbb4','Lamc1','Dag1','Tgfb2','Tgfbr1_Tgfbr2']
 with h5py.File(path) as h:
  obs=read_elem(h['obs']);var=read_elem(h['var']);symbols=var.gene_symbol.astype(str).to_numpy()
  assert obs.index.equals(annotation.index)
  wanted=sorted(set('_'.join(checked+selected.receptor.unique().tolist()).split('_')));ij={g:np.flatnonzero(symbols==g) for g in wanted}
  x=sparse_dataset(h['X'])[:].tocsr();den=np.asarray(x.sum(1)).ravel()
  counts=np.column_stack([np.asarray(x[:,ij[g]].sum(1)).ravel() for g in wanted]);del x
  values=np.log1p(counts*10000/np.maximum(den,1)[:,None])
  expression_checks=0
  for group in ['CM','Broad:Endothelial','Broad:Fibroblast']:
   for stage in STAGES:
    for lib in annotation.loc[annotation.age.eq(stage),'sample'].unique():
     mask=annotation.age.eq(stage)&annotation['sample'].eq(lib)
     mask &= annotation.receiver_ventricular if group=='CM' else annotation.source_broad.eq(group[6:])
     if not mask.any():continue
     for complex_ in checked:
      ids=[wanted.index(g) for g in complex_.split('_')]
      v=values[mask][:,ids].min(axis=1);det=(counts[mask][:,ids]>=1).all(axis=1)
      target=corrected[corrected.stage.eq(stage)&corrected['sample'].eq(lib)&corrected.source_group.eq(group)&corrected.species.eq(complex_)].iloc[0]
      np.testing.assert_allclose(v.mean(),target.mean_log1p_CP10k,rtol=2e-6,atol=2e-8)
      np.testing.assert_allclose(det.mean(),target.fraction_complex_detected,rtol=2e-6,atol=2e-8)
      expression_checks+=1
 specificity=[]
 # The initial non-atrial screen includes provisional conduction candidates.
 # Independently assess whether selected receptors also have support after
 # excluding those candidates and uncharacterized ventricular candidates.
 clean=annotation.receiver_ventricular&~annotation.CM_state.isin(['Conduction-program candidate','Uncharacterized ventricular candidate'])
 for stage in STAGES:
  for lib in annotation.loc[annotation.age.eq(stage),'sample'].unique():
   mask=clean&annotation.age.eq(stage)&annotation['sample'].eq(lib)
   for rec in selected.receptor.unique():
    ids=[wanted.index(g) for g in rec.split('_')]
    mu=values[mask][:,ids].min(axis=1).mean();fraction=(counts[mask][:,ids]>=1).all(axis=1).mean()
    specificity.append({'stage':stage,'sample':lib,'receptor':rec,'non_conduction_CM_nuclei':int(mask.sum()),
     'mean_log1p_CP10k':mu,'fraction_all_subunits_detected':fraction,'library_support':bool(mask.sum()>=30 and mu>=.03 and fraction>=.02)})
 specificity=pd.DataFrame(specificity);save(specificity,'RNA_receptor_nonconduction_CM_sensitivity.tsv')
 qualified=specificity.groupby(['stage','receptor']).library_support.all()
 stagepairs=read(OUT/'RNA_pair_stage_support.tsv');kept=stagepairs[stagepairs.RNA_supported&stagepairs.pair_id.isin(selected.pair_id)].merge(selected[['pair_id','receptor']],on='pair_id')
 kept['nonconduction_receptor_two_library_support']=[bool(qualified.loc[(s,r)]) for s,r in zip(kept.stage,kept.receptor)]
 save(kept,'RNA_selected_pair_receiver_sensitivity.tsv')
 result={'status':'PASS','RNA_support_rows_independently_reconstructed':len(support),
  'RNA_library_complex_means_and_detection_recomputed_from_full_counts':expression_checks,
  'all_selected_subunits_measured':True,'selected_predictor_outcome_gene_overlap':0,
  'non_conduction_CM_receptor_sensitivity_supported_stage_pairs':int(kept.nonconduction_receptor_two_library_support.sum()),
  'non_conduction_CM_receptor_sensitivity_total_stage_pairs':len(kept)}
 (OUT/'RNA_independent_validation.json').write_text(json.dumps(result,indent=2)+'\n');print('JOINT_RNA_VALIDATION',result,flush=True)
 return result

def spatial_audit():
 diagnostic=[];r2checks=0;gains=[];cohort=[];local_checks=0
 for stage in STAGES:
  sources=read(OUT/stage/'source_domains.tsv.gz',index_col=0);receivers=read(OUT/stage/'receiver_domains.tsv.gz',index_col=0)
  assert sources.index.is_unique and receivers.index.is_unique
  full=read(f'results/cardiac_spatial_pipeline/{stage}/atlas_annotations.tsv.gz',index_col=0)
  for fold in range(4):
   train=set(sources.index[sources.domain.ne(fold)]);test=set(sources.index[sources.domain.eq(fold)])
   assert not train&test
  for view in ['core','native']:
   di=read(OUT/stage/(view+'_transport_diagnostics.tsv'))
   assert len(di)==4 and di.converged.all()
   assert (di.maximum_relative_KKT_residual<=1.001e-3).all()
   assert (di.maximum_relative_supply_excess<=1.001e-3).all() and (di.maximum_relative_demand_excess<=1.001e-3).all()
   diagnostic.append(di)
   with np.load(OUT/stage/(view+'_predictor_features.npz')) as z:
    ids=z['cell_ids'];xy=z['xy'];domain=z['domain'];take=z['eligible'];X=z['features'];Y=z['outcomes'];names=z['columns'];yn=z['outcome_names']
   assert np.isfinite(X[take]).all() and np.isfinite(Y[take]).all()
   # Reconstruct local ligand exposure and source-capture controls directly
   # from counts and coordinates, without the production sparse graph.
   matrix=(ROOT/f'results/xenium_nuclear_anchored_segmentation/{stage}/nuclear_anchored_counts.h5ad' if view=='core' else ROOT/f'results/xenium_ventricular_mapping/{stage}/{stage}_xenium_counts_and_candidate_mapping.h5ad')
   with h5py.File(matrix) as h:
    obs=read_elem(h['obs']).index;genes=read_elem(h['var']).index
    xx=sparse_dataset(h['layers/high_confidence'] if view=='core' else h['X'])[obs.get_indexer(sources.index)]
    depth=np.asarray(xx[:,~genes.isin(OUTCOME_GENES)].sum(1)).ravel()
    ligands={g:np.asarray(xx[:,genes.get_loc(g)].toarray()).ravel()*10000/np.maximum(depth,1) for g in ['Lamc1','Nrg1']};del xx
   with np.load(ROOT/f'results/cardiac_spatial_reconstruction/{stage}/morphology_masks.npz') as z:
    tissue=ndimage.binary_dilation(z['tissue'],iterations=1);resolution=float(z['resolution_um'])
   sxy=full.loc[sources.index,['x_centroid','y_centroid']].to_numpy();tree=cKDTree(sxy)
   for fold in range(4):
    samples=np.flatnonzero(take&(domain==fold));samples=samples[np.linspace(0,len(samples)-1,8).astype(int)]
    for k in samples:
     center=full.loc[ids[k],['x_centroid','y_centroid']].to_numpy(float)
     si=np.array(tree.query_ball_point(center,200),int);si=si[(sources.domain.iloc[si].to_numpy()==fold)&(sources.index[si]!=ids[k])]
     keep=np.ones(len(si),bool)
     for weight in [.15,.3,.5,.7,.85]:
      q=np.rint(((1-weight)*sxy[si]+weight*center)[:,::-1]/resolution).astype(int)
      q=np.clip(q,0,np.array(tissue.shape)-1);keep &= tissue[q[:,0],q[:,1]]
     si=si[keep];distance=np.maximum(np.linalg.norm(sxy[si]-center,axis=1),.001)
     own=sources.index.get_loc(ids[k]);np.testing.assert_allclose(X[k,list(names).index('A__log_CM_capture')],np.log1p(depth[own]),atol=2e-6,rtol=2e-6);local_checks+=1
     for source,gene,radius in [('Broad:Fibroblast','Lamc1',50),('Broad:Endothelial','Nrg1',200)]:
      group='FB:' if source.endswith('Fibroblast') else 'EC:'
      chosen=sources.source_group.iloc[si].str.startswith(group).to_numpy()&(distance<=radius)
      ii=si[chosen];weights=np.exp(-distance[chosen]/(radius/2));den=weights.sum()
      control='A__source_capture__'+source+f'__{radius}um'
      expected=(weights*np.log1p(depth[ii])).sum()/max(den,1e-8)
      np.testing.assert_allclose(X[k,list(names).index(control)],expected,atol=5e-6,rtol=5e-6);local_checks+=1
      for stat in ['sum','mean']:
       name=f'ligand__{source}__{gene}__{radius}um__{stat}'
       if name not in names:continue
       expected=(weights*ligands[gene][ii]).sum()
       if stat=='mean':expected/=max(den,1e-8)
       np.testing.assert_allclose(X[k,list(names).index(name)],np.log1p(expected),atol=5e-6,rtol=5e-6);local_checks+=1
   np.testing.assert_array_equal(take,receivers.analysis_eligible.to_numpy())
   frozen=read(f'results/CM_maturation21_nonCM_context/{stage}/CM_maturation21_scores.tsv.gz',index_col='cell_id')
   assert frozen.loc[ids,'strict_ventricular_map_cohort'].all()
   np.testing.assert_allclose(Y,frozen.loc[ids,yn].to_numpy(),rtol=2e-6,atol=1e-6)
   for fold in range(4):
    train=take&(domain!=fold);test=take&(domain==fold)
    gap=cKDTree(xy[test]).query(xy[train])[0].min();assert gap>=399.99
    # Check actual source neighborhoods in the original, unsplit geometry too.
    tt=cKDTree(sources[['x_centroid','y_centroid']].to_numpy())
    for k in np.flatnonzero(test)[::max(1,int(test.sum()/64))]:
     near=tt.query_ball_point(xy[k],199.99)
     assert (sources.domain.iloc[near]==fold).all()
   with np.load(OUT/stage/(view+'_OOF_predictions.npz')) as z:
    P=z['predictions'];sets=list(z['feature_sets']);learners=list(z['learners'])
    np.testing.assert_allclose(Y,z['outcomes']);assert np.isfinite(P[:,:,take]).all();assert np.isnan(P[:,:,~take]).all()
    assert float(z['ridge_minimum_molecular_scale'])==.05
    np.testing.assert_array_equal(z['predictor_columns'],names)
   performance=read(OUT/stage/(view+'_OOF_performance.tsv'))
   for r in performance[performance.CM_state.eq('All')].itertuples():
    j=list(yn).index(r.outcome);li=learners.index(r.learner);si=sets.index(r.feature_set)
    independent=1-np.sum((Y[take,j]-P[li,si,take,j])**2)/np.sum((Y[take,j]-Y[take,j].mean())**2)
    np.testing.assert_allclose(independent,r.OOF_R2,rtol=2e-5,atol=2e-6);r2checks+=1
   for li,learner in enumerate(learners):
    for bigger,smaller in [('D','C_fine_context'),('E','C_fine_context'),('E','D')]:
     for j,key in enumerate(yn):
      delta=r2_score(Y[take,j],P[li,sets.index(bigger),take,j])-r2_score(Y[take,j],P[li,sets.index(smaller),take,j])
      gains.append({'stage':stage,'view':view,'learner':learner,'comparison':bigger+' minus '+smaller,'CM_state':'All','outcome':key,'independent_delta_R2':delta})
   cohort.append({'stage':stage,'view':view,'mapped_CMs':len(ids),'test_eligible_CMs':int(take.sum()),'all_cardiac_source_profiles':len(sources),'features':len(names)})
 gain=pd.DataFrame(gains);reported=read(OUT/'incremental_prediction_gains.tsv')
 merged=gain.merge(reported,on=['stage','view','learner','comparison','CM_state','outcome'],validate='1:1')
 np.testing.assert_allclose(merged.independent_delta_R2,merged.delta_R2,rtol=1e-4,atol=1e-6)
 primary=reported[reported.primary_analysis&reported.blocks800.ge(8)].sort_values('within_section_p')
 p=primary.within_section_p.to_numpy();q=np.minimum(1,np.minimum.accumulate((p*len(p)/np.arange(1,len(p)+1))[::-1])[::-1])
 np.testing.assert_allclose(q,primary.primary_grid_BH_q,rtol=2e-5,atol=1e-8)
 save(pd.concat(diagnostic),'all_transport_diagnostics.tsv');save(pd.DataFrame(cohort),'analysis_cohorts.tsv')
 result={'status':'PASS','transport_domains_checked':sum(len(z) for z in diagnostic),'OOF_R2_independently_recomputed':r2checks,
  'incremental_prediction_gains_recomputed':len(merged),'primary_BH_tests':len(primary),'source_ID_overlap_across_domains':0,
  'frozen_observed_CM_outcomes_verified':True,'minimum_receiver_gap_um':400,'local_source_sets_independently_checked':True,
  'all_completed_transport_plans_meet_KKT_tolerance':True,'new_independent_spatial_hearts':0}
 result['independently_recomputed_local_ligand_and_capture_checks']=local_checks
 (OUT/'spatial_independent_validation.json').write_text(json.dumps(result,indent=2)+'\n');print('JOINT_SPATIAL_VALIDATION',result,flush=True)
 return result

if __name__=='__main__':
 import sys
 rna_audit() if '--RNA-only' in sys.argv else spatial_audit()
 export_inputs()
