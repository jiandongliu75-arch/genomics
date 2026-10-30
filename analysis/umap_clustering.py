"""Cluster an existing expression-neighbor graph and export GitHub-sized analysis results."""
from pathlib import Path
import os,argparse,json,hashlib,importlib.metadata,tempfile
parser=argparse.ArgumentParser()
parser.add_argument('--input',type=Path,required=True)
parser.add_argument('--annotations',type=Path)
parser.add_argument('--outdir',type=Path,required=True)
args=parser.parse_args();out=args.outdir;out.mkdir(parents=True,exist_ok=True)
os.environ.setdefault('MPLCONFIGDIR',str(Path(tempfile.gettempdir())/'genomics-umap-matplotlib'))
import numpy as np,pandas as pd,scanpy as sc,anndata as ad
from scipy import sparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sc.settings.n_jobs=4
print('READ',flush=True)
a=sc.read_h5ad(args.input,backed='r');obs=a.obs.copy();var=a.var.copy();xy=a.obsm['X_umap'].copy()
q=ad.AnnData(sparse.csr_matrix((a.n_obs,0)),obs=obs.copy())
q.obsp['connectivities']=a.obsp['connectivities'].copy()
q.obsp['distances']=a.obsp['distances'].copy();q.uns['neighbors']=dict(a.uns['neighbors'])
q.obsm['X_umap']=xy;q.obsm['X_pca']=a.obsm['X_pca'].copy()
for resolution in [.4,.8,1.2]:
 key=f'leiden_r{resolution:g}'
 sc.tl.leiden(q,resolution=resolution,random_state=42,flavor='igraph',n_iterations=-1,key_added=key)
 print(key,q.obs[key].nunique(),flush=True)
q.obs['cluster']=q.obs['leiden_r0.8'].copy()
if args.annotations:
 ann=pd.read_csv(args.annotations,index_col=0)
 assert ann.index.is_unique and set(q.obs_names)==set(ann.index)
 q.obs['preliminary_lineage']=ann.loc[q.obs_names,'preliminary_lineage'].astype(str)
else:q.obs['preliminary_lineage']='Not annotated'
# Publish actual coordinates and all QC metadata, not merely rendered images.
q.obs['UMAP1']=xy[:,0];q.obs['UMAP2']=xy[:,1]
q.obs.to_csv(out/'cell_clusters_and_umap.csv.gz',compression={'method':'gzip','compresslevel':1})
clusters=sorted(q.obs.cluster.unique(),key=int);cluster_index=q.obs.cluster.map({c:i for i,c in enumerate(clusters)}).to_numpy().astype(int)
group=sparse.csr_matrix((np.ones(a.n_obs,dtype=np.float64),(cluster_index,np.arange(a.n_obs))),shape=(len(clusters),a.n_obs))
n=np.asarray(group.sum(axis=1)).ravel();assert (n>0).all()
print('EXPRESSION_SUMMARIES',flush=True)
x=a.X[:].tocsr().astype(np.float32)
# Fractional pseudobulk values preserve the DecontX output without rounding.
sums=(group@x).toarray()
expr=pd.DataFrame(sums.T,index=var.index,columns=clusters)
expr.index.name='gene_id';expr.insert(0,'gene_symbol',var.gene_symbol.to_numpy())
expr.to_csv(out/'cluster_pseudobulk_corrected_counts.tsv.gz',sep='\t',float_format='%.6g',compression={'method':'gzip','compresslevel':1})
# Export lineage marker counts per nucleus for plotting/annotation review.
markers=['Tnnt2','Tnni3','Actc1','Myh6','Myh7','Ttn','Ryr2','Pdgfra','Dcn','Lum','Col1a1','Col1a2','Col3a1','Cdh5','Pecam1','Kdr','Flt1','Esam','Emcn','Pdgfrb','Rgs5','Cspg4','Kcnj8','Abcc9','Myh11','Acta2','Tagln','Cnn1','Ptprc','Tyrobp','Lyz2','Csf1r','Cd3d','Cd79a','Upk3b','Wt1','Msln','Krt19','Plp1','Sox10','Mpz','S100b']
idx=np.flatnonzero(var.gene_symbol.isin(markers))
marker_frame=pd.DataFrame(x[:,idx].toarray(),index=q.obs_names,columns=var.iloc[idx].gene_symbol.to_numpy())
marker_frame.index.name='cell_id';marker_frame.to_csv(out/'lineage_marker_corrected_counts.csv.gz',float_format='%.5g',compression={'method':'gzip','compresslevel':1})
# Descriptive marker ranking on log-normalized expression; no p-values implied.
detected=x.copy();detected.data[:]=1
fraction=(group@detected).toarray()/n[:,None]
del detected
x=x.multiply(1e4/q.obs.total_counts.to_numpy()[:,None]).tocsr();x.data=np.log1p(x.data)
logsum=(group@x).toarray();total_log=np.asarray(x.sum(axis=0)).ravel();total_detect=(fraction*n[:,None]).sum(axis=0)
ranked=[];summary=[]
for i,c in enumerate(clusters):
 mean=logsum[i]/n[i];rest=(total_log-logsum[i])/(a.n_obs-n[i]);diff=mean-rest
 frest=(total_detect-fraction[i]*n[i])/(a.n_obs-n[i])
 eligible=(fraction[i]>=.2)&((fraction[i]-frest)>=.1)&(~var.mt.to_numpy())&(~var.ribo.to_numpy())
 order=np.flatnonzero(eligible);order=order[np.argsort(diff[order])[::-1]][:30]
 for rank,j in enumerate(order,1):ranked.append({'cluster':c,'rank':rank,'gene_id':var.index[j],'gene_symbol':var.gene_symbol.iloc[j],'mean_log1p_10k_difference':float(diff[j]),'fraction_detected_cluster':float(fraction[i,j]),'fraction_detected_rest':float(frest[j])})
 sub=q.obs[q.obs.cluster==c];labels=sub.preliminary_lineage.value_counts();label=labels.index[0];purity=labels.iloc[0]/len(sub)
 summary.append({'cluster':c,'nuclei':len(sub),'dominant_preliminary_lineage':label,'lineage_fraction':purity,'median_estimated_contamination':sub.decontx_contamination.median(),'median_corrected_umis':sub.total_counts.median(),'top_descriptive_markers':';'.join(var.gene_symbol.iloc[order[:10]]),'small_cluster_flag':len(sub)<50,'mixed_lineage_flag':purity<.8})
pd.DataFrame(ranked).to_csv(out/'cluster_top_markers.tsv',sep='\t',index=False)
sumframe=pd.DataFrame(summary);sumframe.to_csv(out/'cluster_summary.tsv',sep='\t',index=False)
for col in ['age','sample','preliminary_lineage']:
 pd.crosstab(q.obs.cluster,q.obs[col]).reindex(clusters).to_csv(out/f'cluster_counts_by_{col}.tsv',sep='\t')
# A metadata-only H5AD is intentionally small enough for ordinary GitHub files.
q.uns['clustering_method']={'algorithm':'Leiden on expression-space kNN graph; not UMAP coordinates','primary_resolution':.8,'sensitivity_resolutions':[.4,.8,1.2],'seed':42,'iterations':'until convergence','batch_integration':'none','expression_matrix':'Absent in this metadata-only object; use counts tables or original full H5AD','lineage_labels':'Earlier preliminary broad labels; not final biological identities'}
q.write_h5ad(out/'umap_clusters_metadata_only.h5ad',compression='gzip')
cmap=plt.get_cmap('tab20');palette={c:cmap(i%20) for i,c in enumerate(clusters)}
fig,axes=plt.subplots(1,2,figsize=(17,7))
for c in clusters:
 mask=q.obs.cluster.eq(c).to_numpy()
 for ax in axes:ax.scatter(xy[mask,0],xy[mask,1],s=1,color=palette[c],rasterized=True)
 center=np.median(xy[mask],axis=0)
 axes[0].text(*center,c,fontsize=9,ha='center',va='center',bbox={'facecolor':'white','alpha':.8,'pad':1,'edgecolor':'none'})
axes[0].set_title(f'Leiden clusters (resolution 0.8; {len(clusters)} clusters)')
axes[1].clear()
for age in ['P0','P7','P14','P21']:
 mask=q.obs.age.eq(age).to_numpy();axes[1].scatter(xy[mask,0],xy[mask,1],s=1,label=age,rasterized=True)
axes[1].legend(markerscale=6);axes[1].set_title('Same nuclei colored by age (no integration)')
for ax in axes:ax.set(xlabel='UMAP1',ylabel='UMAP2')
fig.tight_layout();fig.savefig(out/'umap_leiden_clusters.png',dpi=180);fig.savefig(out/'umap_leiden_clusters.pdf');plt.close(fig)
fig,axes=plt.subplots(1,3,figsize=(18,5))
for ax,res in zip(axes,[.4,.8,1.2]):
 col=f'leiden_r{res:g}'
 for j,c in enumerate(sorted(q.obs[col].unique(),key=int)):
  mask=q.obs[col].eq(c).to_numpy();ax.scatter(xy[mask,0],xy[mask,1],s=.7,color=cmap(j%20),rasterized=True)
 ax.set(title=f'Resolution {res}: {q.obs[col].nunique()} clusters',xlabel='UMAP1',ylabel='UMAP2')
fig.tight_layout();fig.savefig(out/'clustering_resolution_comparison.png',dpi=160);fig.savefig(out/'clustering_resolution_comparison.pdf');plt.close(fig)
# Validate round-trip IDs and the published data dimensions.
check=pd.read_csv(out/'cell_clusters_and_umap.csv.gz',index_col=0)
assert len(check)==a.n_obs and check.index.is_unique and check.cluster.notna().all()
assert np.isfinite(check[['UMAP1','UMAP2']]).all().all()
assert int(sumframe.nuclei.sum())==a.n_obs
assert expr.shape[0]==a.n_vars and np.allclose(sums.sum(axis=0),np.asarray(a.X[:].astype(np.float64).sum(axis=0)).ravel(),rtol=1e-5,atol=.05)
a.file.close()
versions={p:importlib.metadata.version(p) for p in ['scanpy','anndata','numpy','pandas','scipy','matplotlib','igraph','leidenalg']}
(out/'requirements.txt').write_text('\n'.join(k+'=='+v for k,v in versions.items())+'\n')
def sha(path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  for chunk in iter(lambda:f.read(1048576),b''):h.update(chunk)
 return h.hexdigest()
manifest={'input_h5ad_sha256':sha(args.input),'nuclei':len(q.obs),'genes_in_input':len(var),'clusters':len(clusters),'algorithm':dict(q.uns['clustering_method']),'versions':versions,'files':[]}
for p in sorted(out.iterdir()):
 if p.is_file() and p.name!='manifest.json':manifest['files'].append({'file':p.name,'bytes':p.stat().st_size,'sha256':sha(p)})
assert all(p['bytes']<90_000_000 for p in manifest['files'])
(out/'manifest.json').write_text(json.dumps(manifest,indent=2))
print('COMPLETE',len(q.obs),'nuclei',len(clusters),'clusters',flush=True)
print(sumframe[['cluster','nuclei','dominant_preliminary_lineage','lineage_fraction','small_cluster_flag']].to_string(index=False),flush=True)
