"""Check corrected-count state markers within age and separately in its two libraries."""
from pathlib import Path
import argparse
import numpy as np,pandas as pd,h5py
from anndata.io import read_elem,sparse_dataset
from scipy import sparse
p=argparse.ArgumentParser();p.add_argument('--input',default='/workspace/heart-atlas/results/cell_identities_subtypes/snrna_qc_annotated_subclustered.h5ad');args=p.parse_args()
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/ventricular_stage_batch'
ann=pd.read_csv(OUT/'ventricular_cells_embeddings_and_clusters.csv.gz',index_col=0)
with h5py.File(args.input) as f:
 obs=read_elem(f['obs']);var=read_elem(f['var']);ix=obs.index.get_indexer(ann.index);assert (ix>=0).all();x=sparse_dataset(f['X'])[ix,:].tocsr().astype(np.float32)
log=(sparse.diags(1e4/np.asarray(x.sum(axis=1)).ravel())@x).tocsr();log.data=np.log1p(log.data)
def mean(m):return np.asarray(log[m].mean(axis=0)).ravel()
def frac(m):return np.asarray((x[m]>=1).mean(axis=0)).ravel()
sexgenes=['Xist','Tsix','Ddx3y','Kdm5d','Eif2s3y','Uty','Usp9y','Zfy1','Zfy2','Sry','Rbmy']
valid=~var.mt.to_numpy()&~var.ribo.to_numpy()&~var.gene_symbol.isin(sexgenes).to_numpy()
rows=[];checks=[]
for cluster,s in ann.groupby('ventricular_cluster'):
 stage=s.age.value_counts().index[0];base=ann.age.eq(stage).to_numpy();a=base&ann.ventricular_cluster.eq(cluster).to_numpy();b=base&~a
 if min(a.sum(),b.sum())<50:continue
 delta=mean(a)-mean(b);fa=frac(a);fb=frac(b);libs=[];contrasts=[]
 for lib in ann.loc[base,'sample'].unique():
  sa=a&ann['sample'].eq(lib).to_numpy();sb=b&ann['sample'].eq(lib).to_numpy()
  if min(sa.sum(),sb.sum())>=20:libs.append(lib);contrasts.append(mean(sa)-mean(sb))
 supported=np.sum(np.array(contrasts)>.05,axis=0) if contrasts else np.zeros(x.shape[1],dtype=int)
 eligible=valid&(delta>.1)&(fa>=.2)&((fa-fb)>=.05);order=np.flatnonzero(eligible);order=order[np.argsort(delta[order])[::-1]][:40]
 for rank,j in enumerate(order,1):
  row={'ventricular_cluster':cluster,'dominant_stage':stage,'rank':rank,'gene_id':var.index[j],'gene_symbol':var.gene_symbol.iloc[j],'mean_log1p_difference_vs_same_stage_rest':delta[j],'fraction_corrected_ge1_cluster':fa[j],'fraction_corrected_ge1_stage_rest':fb[j],'eligible_libraries':len(libs),'libraries_with_difference_gt_0.05':int(supported[j])}
  for lib,d in zip(libs,contrasts):row[lib+'_difference']=d[j]
  rows.append(row)
 supportedgenes=[str(var.gene_symbol.iloc[j]) for j in order if supported[j]>=2]
 checks.append({'ventricular_cluster':cluster,'dominant_stage':stage,'n_in_dominant_stage':int(a.sum()),'eligible_libraries':len(libs),'top40_markers_supported_in_both_libraries':len(supportedgenes),'supported_marker_names':';'.join(supportedgenes[:12])})
pd.DataFrame(rows).to_csv(OUT/'within_stage_descriptive_markers.tsv',sep='\t',index=False)
pd.DataFrame(checks).to_csv(OUT/'within_stage_marker_library_support.tsv',sep='\t',index=False)
print(pd.DataFrame(checks).to_string(index=False))
