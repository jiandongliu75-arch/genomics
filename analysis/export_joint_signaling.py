"""Compact per-cell source data and descriptive molecular availability tables."""
import json
import numpy as np
import pandas as pd
from joint_signaling_common import *

MAP_FEATURES=[
 'ligand__Broad:Endothelial__Nrg1__200um__sum','ligand__Broad:Fibroblast__Lamc1__50um__sum',
 'cot__Broad:Fibroblast__Lamc1--Dag1','cot__Broad:Fibroblast__Fn1--Itga5_Itgb1',
 'cot__Broad:Fibroblast__Tnc--Itga9_Itgb1','cot__Broad:Endothelial__Nrg1--Erbb4',
 'cot__Broad:Endothelial__Lama4--Dag1','cot__Broad:Endothelial__Igf2--Igf1r',
 'cot__Broad:Endothelial__Ptprm--Ptprm','cot__CM__Ptprm--Ptprm',
 'cot__FB:Postn__Tnc--Itga9_Itgb1','cot__FB:Col8a1__Lamc1--Dag1',
 'cot__FB:Gfpt2__Lamc1--Dag1','cot__FB:Lsamp__Lamc1--Dag1',
 'cot__EC:Endocardial-like__Nrg1--Erbb4','cot__EC:Capillary-like__Nrg1--Erbb4',
 'receptor__Erbb4','receptor__Dag1','receptor__Itga5_Itgb1']

def main():
 programs=read(OUT/'RNA_CM_programs_by_library.tsv')
 programs=programs[programs.CM_state.isin(['Ankrd1-associated','Other working CM'])]
 contrast=programs.pivot(index=['stage','sample','program'],columns='CM_state',values='mean_log1p_CP10k').reset_index()
 contrast['Ankrd1_minus_working']=contrast['Ankrd1-associated']-contrast['Other working CM']
 save(contrast,'RNA_Ankrd1_working_independent_program_contrasts.tsv')
 expr=read(OUT/'RNA_complex_expression_by_library.tsv.gz');expr=expr[expr.view.eq('corrected')]
 support=read(OUT/'RNA_pair_source_support.tsv.gz')
 preferred=[('Broad:Endothelial','Nrg1--Erbb4'),('Broad:Endothelial','Igf2--Igf1r'),
 ('Broad:Fibroblast','Fn1--Itga5_Itgb1'),('Broad:Fibroblast','Lamc1--Dag1'),
 ('Broad:Fibroblast','Tnc--Itga9_Itgb1'),('CM','Ptprm--Ptprm')]
 rows=[]
 for group,pair in preferred:
  for stage in STAGES:
   r=support[support.source_group.eq(group)&support.pair_id.eq(pair)&support.stage.eq(stage)].iloc[0]
   src=expr[expr.source_group.eq(group)&expr.stage.eq(stage)&expr.species.eq(r.ligand)]
   rec=expr[expr.source_group.eq('CM')&expr.stage.eq(stage)&expr.species.eq(r.receptor)]
   rows.append({'stage':stage,'source':group,'pair_id':pair,'RNA_two_library_availability':r.joint_two_library_support,
    'source_mean_log1p_CP10k':src.mean_log1p_CP10k.mean(),'receiver_mean_log1p_CP10k':rec.mean_log1p_CP10k.mean(),
    'source_detection_min':src.fraction_complex_detected.min(),'source_detection_max':src.fraction_complex_detected.max(),
    'receiver_detection_min':rec.fraction_complex_detected.min(),'receiver_detection_max':rec.fraction_complex_detected.max(),
    'evidence':'Corrected nuclear RNA availability; not signaling activation or source-specific causality'})
 save(pd.DataFrame(rows),'biologically_selected_RNA_availability.tsv')
 summaries=[]
 for stage in STAGES:
  for view in ['core','native']:
   with np.load(OUT/stage/(view+'_predictor_features.npz')) as z:
    ids=z['cell_ids'];take=z['eligible'];columns=list(z['columns']);values=z['features'];xy=z['xy'];Y=z['outcomes'];yn=z['outcome_names']
    d=pd.DataFrame({'cell_id':ids,'stage':stage,'view':view,'CM_state':z['CM_states'],'x_centroid':xy[:,0],
       'y_centroid':xy[:,1],'domain':z['domain'],'eligible_spatial_test':take})
   for j,name in enumerate(yn):d[name]=Y[:,j]
   for feature in MAP_FEATURES:
    d[feature]=values[:,columns.index(feature)] if feature in columns else np.nan
    if feature in columns:
     v=values[take,columns.index(feature)]
     summaries.append({'stage':stage,'view':view,'feature':feature,'test_eligible_CMs':int(take.sum()),
       'mean_log_input':v.mean(),'median_log_input':np.median(v),'fraction_positive':(v>1e-12).mean(),
       'interpretation':'Observed receptor log-expression or log1p model-allocated candidate input; not activated signaling'})
   with np.load(OUT/stage/(view+'_OOF_predictions.npz')) as z:
    P=z['predictions'];sets=list(z['feature_sets']);learners=list(z['learners'])
    for feature_set in ['A','C_fine_context','D','E']:
     for j,name in enumerate(yn):d['OOF_RF_'+feature_set+'__'+name]=P[learners.index('RandomForest'),sets.index(feature_set),:,j]
   save(d,stage+'/'+view+'_CM_candidate_map_values.tsv.gz')
 save(pd.DataFrame(summaries),'selected_spatial_input_summary.tsv')
 gains=read(OUT/'incremental_prediction_gains.tsv');primary=gains[gains.primary_analysis]
 perf=read(OUT/'OOF_performance.tsv')
 q=primary[primary.primary_grid_BH_q.lt(.05)]
 summary={'primary_supported_tests':int(primary.primary_grid_BH_q.notna().sum()),'primary_BH_significant_positive_gains':int(q.delta_R2.gt(0).sum()),
  'primary_BH_significant_negative_gains':int(q.delta_R2.lt(0).sum()),'primary_gain_min':primary.delta_R2.min(),'primary_gain_max':primary.delta_R2.max(),
  'interpretation':'Prediction gains concern fixed models in unseen regions of the same section, not causal effects or independent animals'}
 (OUT/'result_summary.json').write_text(json.dumps(summary,indent=2)+'\n');export_inputs();print('JOINT_SOURCE_EXPORTS',summary,flush=True)
if __name__=='__main__':main()
