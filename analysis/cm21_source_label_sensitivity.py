"""Challenge source annotations while holding the 21-gene CM outcome fixed."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','2')
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.spatial import cKDTree
from cm21_paperstyle_extension import ROOT,OUT,STAGES,Design
import cm21_primary_analysis as old
def read(p,**kw):return pd.read_csv(ROOT/p,sep='\t',**kw)
def main():
 rows=[];coverage=[]
 for stage in STAGES:
  d=read(f'results/CM21_primary_spatial/{stage}_CM21_scores_and_neighborhoods.tsv.gz',index_col=0)
  extra=read(f'results/CM21_paperstyle_extension/{stage}_CM21_extension_values.tsv.gz',index_col=0)
  for col in ['native_module21','core_module21','native_excluded21_depth','core_excluded21_depth']:d[col]=extra[col]
  f=read(f'results/cardiac_spatial_pipeline/{stage}/atlas_annotations.tsv.gz',index_col=0)
  base,catalog=old.source_masks(stage,f);available=set(catalog[catalog.support_status.eq('Available')].source)
  ec=f.atlas_identity.eq('Endothelial')&f.program_source_eligible
  masks={('EC:Capillary-like','reviewed_nuclear_supported'):base['EC:Capillary-like'],
   ('EC:Capillary-like','public_capillary_with_EC_QC'):ec&f.author_atlas_cell_state.eq('cap.EC'),
   ('EC:Capillary-like','strict_cell_and_nuclear'):ec&f.fine_identity.eq('Capillary-like')}
  fb=read(f'results/fibroblast_multimethod_identity/{stage}_multimethod_FB_annotations.tsv.gz',index_col=0)
  for family in old.FB:
   source='FB:'+family
   if source not in available:continue
   for layer in ['recommended_candidate','crosschecked_candidate','boundary_supported_candidate','cell_supported_candidate']:
    mask=pd.Series(False,index=f.index);mask.loc[fb.index]=fb.program_source_eligible&fb[layer].eq(family+'-associated')
    masks[(source,layer)]=mask&base['Fibroblast']
  xy=f[['x_centroid','y_centroid']].to_numpy();ci=f.index.get_indexer(d.index)
  G=sparse.load_npz(ROOT/f'results/cardiac_spatial_pipeline/{stage}/proximity50_graph.npz').tocsr()[ci]
  designs={v:Design(d,v) for v in ['native','core']}
  for (source,layer),mask in masks.items():
   assert not d.index.isin(f.index[mask]).any()
   count=int(mask.sum())
   if count<30:
    coverage.append(dict(stage=stage,source=source,layer=layer,sources=count,assessable=False));continue
   distance=cKDTree(xy[mask]).query(xy[ci])[0];near=distance<=100
   blocks=len(np.unique(xy[ci][near]//400,axis=0));supported=int(near.sum())>=100 and blocks>=8
   coverage.append(dict(stage=stage,source=source,layer=layer,sources=count,CMs_within100=int(near.sum()),nearby_blocks=blocks,assessable=supported))
   if not supported:continue
   xx={'nearness':-np.log1p(distance)}
   if source.startswith('EC:'):xx['count50']=np.log1p(np.asarray(G@mask.to_numpy(float)).ravel())
   for view,design in designs.items():
    for metric,x in xx.items():
     rr=design.fit(x,d[[view+'_module21']].to_numpy(),['module_ensemble'])
     rows.extend(dict(stage=stage,view=view,source=source,layer=layer,metric=metric,sources=count,**r) for r in rr)
 pd.DataFrame(rows).to_csv(OUT/'source_label_sensitivity.tsv',sep='\t',index=False,float_format='%.9g')
 pd.DataFrame(coverage).to_csv(OUT/'source_label_sensitivity_coverage.tsv',sep='\t',index=False)
 print(pd.DataFrame(rows).query("source=='EC:Capillary-like' and view=='native' and metric=='nearness'").to_string(index=False))
if __name__=='__main__':main()
