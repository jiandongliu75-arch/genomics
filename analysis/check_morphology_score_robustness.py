"""Targeted score-definition sensitivities for two named follow-up questions."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','2')
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.sparse.csgraph import dijkstra
from scipy.stats import false_discovery_control
from morphology_nonCM_associations import ROOT,OUT,fit,zscore
def main():
 rows=[]
 for stage,source in [('P0','FB:Postn'),('P7','EC:Capillary-like')]:
  d=pd.read_csv(OUT/f'{stage}_CM_region_annotations_and_scores.tsv.gz',sep='\t',index_col=0)
  f=pd.read_csv(ROOT/f'results/cardiac_spatial_pipeline/{stage}/atlas_annotations.tsv.gz',sep='\t',index_col=0)
  old=pd.read_csv(ROOT/f'results/CM_maturation21_nonCM_context/{stage}/CM_maturation21_scores.tsv.gz',sep='\t',index_col=0)
  if stage=='P0':
   fb=pd.read_csv(ROOT/f'results/fibroblast_multimethod_identity/{stage}_multimethod_FB_annotations.tsv.gz',sep='\t',index_col=0)
   mask=pd.Series(False,index=f.index);mask.loc[fb.index]=fb.program_source_eligible&fb.recommended_candidate.eq('Postn-associated')
  else:mask=f.atlas_identity.eq('Endothelial')&f.program_source_eligible&f.EC_nuclear_type.eq('Capillary-like')
  xy=f[['x_centroid','y_centroid']].to_numpy();g=sparse.load_npz(ROOT/f'results/cardiac_spatial_reconstruction/{stage}/tissue_screened_Delaunay_graph.npz').tocoo()
  weighted=sparse.csr_matrix((np.linalg.norm(xy[g.row]-xy[g.col],axis=1),(g.row,g.col)),shape=g.shape)
  path=dijkstra(weighted,directed=False,indices=np.flatnonzero(mask),min_only=True)[f.index.get_indexer(d.index)]
  d['E']=zscore(-np.log1p(path));d['source_local_support']=path<=100
  for name,col in [('no_Myl2','core_Z__without_Myl2'),('without_three_dominant_genes','core_Z__without_Cebpb_Rxrg_Lpin1'),('equal_50_molecules','core_Z__depth50_mean3')]:
   d['S__'+name]=zscore(old.loc[d.index,col])
  target=['Y__core_maturation21','S__no_Myl2','S__without_three_dominant_genes']
  rr,_=fit(d[d.annotated_primary],'E',target,'Score_omissions')
  rows += [dict(stage=stage,source=source,**r) for r in rr]
  deep=d[d.annotated_primary&np.isfinite(d.S__equal_50_molecules)]
  rr,_=fit(deep,'E',['Y__core_maturation21','S__equal_50_molecules'],'Same_depth_eligible_CMs')
  rows += [dict(stage=stage,source=source,**r) for r in rr]
 result=pd.DataFrame(rows);result['score_sensitivity_BH_q']=false_discovery_control(result.within_section_p)
 result.to_csv(OUT/'targeted_score_definition_sensitivities.tsv',sep='\t',index=False,float_format='%.9g')
 print(result[['stage','source','scope','outcome','beta_SD','CI_low','CI_high','score_sensitivity_BH_q']].to_string(index=False))
if __name__=='__main__':main()
