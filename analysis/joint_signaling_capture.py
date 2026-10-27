"""View-specific capture controls, independent of the modeled outcome genes."""
import h5py
import numpy as np
import pandas as pd
from anndata.io import sparse_dataset,read_elem
from joint_signaling_common import *
from joint_signaling_spatial import graph

def augment(stage,view):
 path=OUT/stage/(view+'_predictor_features.npz')
 with np.load(path) as z:
  if 'capture_controls_version' in z and int(z['capture_controls_version'])==1:return
  content={k:z[k] for k in z.files}
 sources=read(OUT/stage/'source_domains.tsv.gz',index_col=0)
 matrix=(ROOT/f'results/xenium_nuclear_anchored_segmentation/{stage}/nuclear_anchored_counts.h5ad' if view=='core' else ROOT/f'results/xenium_ventricular_mapping/{stage}/{stage}_xenium_counts_and_candidate_mapping.h5ad')
 with h5py.File(matrix) as h:
  ids=read_elem(h['obs']).index;genes=read_elem(h['var']).index;ij=ids.get_indexer(sources.index);assert (ij>=0).all()
  x=sparse_dataset(h['layers/high_confidence'] if view=='core' else h['X'])[ij]
  depth=np.asarray(x[:,~genes.isin(OUTCOME_GENES)].sum(1)).ravel();del x
 X=pd.DataFrame(content['features'].astype(float),columns=content['columns'],index=content['cell_ids'])
 ci=sources.index.get_indexer(X.index);valid=ci>=0
 X.loc[valid,'A__log_CM_capture']=np.log1p(depth[ci[valid]])
 names=GROUPS+['Broad:Endothelial','Broad:Fibroblast']
 extra=[]
 for source in names:
  isbroad=not source.startswith(('EC:','FB:'))
  for radius in [50,200]:
   name=('A__source_capture__' if isbroad else 'capture__')+source+f'__{radius}um'
   X[name]=0.
   extra.append({'name':name,'feature_group':'A' if isbroad else 'C_fine_context','source':source,'radius_um':radius,
     'stat':'Kernel-weighted mean log1p outcome-excluded source depth'})
 with np.load(ROOT/f'results/cardiac_spatial_reconstruction/{stage}/morphology_masks.npz') as z:
  tissue=z['tissue'];res=float(z['resolution_um'])
 for domain in range(4):
  si=np.flatnonzero(sources.domain.eq(domain));local=sources.iloc[si]
  M,_=graph(local[['x_centroid','y_centroid']].to_numpy(),tissue,res)
  receiver=X.index.get_indexer(local.index);ok=receiver>=0
  for radius in [50,200]:
   K=M.copy();K.data=np.where(K.data<=radius,np.exp(-K.data/(radius/2)),0);K.eliminate_zeros()
   for source in names:
    mask=(local.source_group.str.startswith('EC:') if source=='Broad:Endothelial' else local.source_group.str.startswith('FB:') if source=='Broad:Fibroblast' else local.source_group.eq(source)).to_numpy(float)
    den=np.asarray(K.T@mask).ravel();num=np.asarray(K.T@(mask*np.log1p(depth[si]))).ravel()
    name=('A__source_capture__' if not source.startswith(('EC:','FB:')) else 'capture__')+source+f'__{radius}um'
    X.loc[X.index[receiver[ok]],name]=(num/np.maximum(den,1e-8))[ok]
 definitions=read(OUT/stage/'feature_definitions.tsv');definitions=definitions[~definitions.name.isin([r['name'] for r in extra])]
 save(pd.concat([definitions,pd.DataFrame(extra)],ignore_index=True),stage+'/feature_definitions.tsv')
 content.update(features=X.to_numpy(np.float32),columns=X.columns.to_numpy(str),capture_controls_version=np.array(1))
 np.savez_compressed(path,**content)
 print('JOINT_CAPTURE_CONTROLS',stage,view,len(extra),flush=True)
