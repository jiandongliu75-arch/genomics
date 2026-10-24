"""Ventricular-only embeddings with stage-protected within-stage sample adjustment."""
from pathlib import Path
import os,argparse,json,hashlib,importlib.metadata
os.environ.setdefault('MPLCONFIGDIR','/tmp/heart-mpl');os.environ.setdefault('NUMBA_CACHE_DIR','/tmp/heart-numba')
import numpy as np,pandas as pd,h5py,scanpy as sc,anndata as ad,harmonypy
from anndata.io import read_elem,sparse_dataset
from scipy import sparse
from sklearn.neighbors import NearestNeighbors
from sklearn.metrics import adjusted_rand_score
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
p=argparse.ArgumentParser();p.add_argument('--input',default='/workspace/heart-atlas/results/cell_identities_subtypes/snrna_qc_annotated_subclustered.h5ad');args=p.parse_args()
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/ventricular_stage_batch';OUT.mkdir(exist_ok=True);sc.settings.n_jobs=4
with h5py.File(args.input) as f:
 allobs=read_elem(f['obs']);var=read_elem(f['var']);indices=np.flatnonzero(allobs.proposed_subtype.eq('Ventricular-like').to_numpy());obs=allobs.iloc[indices].copy()
 print('READ',len(obs),flush=True);x=sparse_dataset(f['X'])[indices,:].tocsr().astype(np.float32)
assert len(obs)==43249 and set(obs.age)=={'P0','P7','P14','P21'}
qc=obs.groupby(['age','sample'],observed=True).agg(nuclei=('age','size'),median_raw_counts=('raw_total_counts','median'),median_raw_genes=('raw_n_genes_by_counts','median'),median_contamination=('decontx_contamination','median'),median_mito=('raw_pct_counts_mt','median'))
qc.to_csv(OUT/'sample_qc_summary.tsv',sep='\t')
work=ad.AnnData(x.copy(),obs=obs[['age','sample']].copy(),var=var.copy());sc.pp.normalize_total(work,target_sum=1e4);sc.pp.log1p(work)
print('HVG_PCA',flush=True);sc.pp.highly_variable_genes(work,n_top_genes=2500,flavor='seurat',batch_key='sample')
sexgenes=['Xist','Tsix','Ddx3y','Kdm5d','Eif2s3y','Uty','Usp9y','Zfy1','Zfy2','Sry','Rbmy']
hv=work.var.highly_variable.to_numpy()&~var.mt.to_numpy()&~var.ribo.to_numpy()&~var.gene_symbol.isin(sexgenes).to_numpy()
hvg=work.var.loc[hv].copy();hvg.to_csv(OUT/'embedding_genes.tsv',sep='\t')
tmp=work[:,hv].copy();sc.pp.scale(tmp,zero_center=False,max_value=10);sc.pp.pca(tmp,n_comps=30,random_state=42,svd_solver='arpack');original=tmp.obsm['X_pca'].copy()
# Run Harmony independently within each stage, in the common global PCA basis.
# Restore each original stage centroid; do not align ages with each other.
adjusted=original.copy();stage_means=[]
for stage in ['P0','P7','P14','P21']:
 m=obs.age.eq(stage).to_numpy();z=original[m];print('HARMONY',stage,len(z),flush=True)
 h=harmonypy.run_harmony(z,obs.loc[m,['sample']],['sample'],theta=2,lamb=1,nclust=30,max_iter_harmony=20,max_iter_kmeans=10,random_state=42,ncores=4,verbose=False)
 corr=np.asarray(h.Z_corr);corr=corr if corr.shape==z.shape else corr.T
 assert corr.shape==z.shape
 corr=corr-corr.mean(axis=0)+z.mean(axis=0);adjusted[m]=corr
 stage_means.append({'stage':stage,'max_centroid_difference':float(np.max(abs(corr.mean(axis=0)-z.mean(axis=0))))})
pd.DataFrame(stage_means).to_csv(OUT/'stage_centroid_validation.tsv',sep='\t',index=False)
# Stronger sensitivity only: within-stage linear removal of measured QC/sex-marker covariates.
# This can remove biological differences and is not the primary representation.
genes=var.gene_symbol.astype(str).to_numpy();sexmarker=[]
for gene in ['Xist','Ddx3y','Kdm5d','Eif2s3y','Uty']:
 j=np.flatnonzero(genes==gene);sexmarker.append(np.asarray(work.X[:,j].sum(axis=1)).ravel())
sexmarker=np.array(sexmarker).T
cov=np.column_stack([np.log1p(obs.raw_total_counts.to_numpy()),np.log1p(obs.raw_n_genes_by_counts.to_numpy()),obs.decontx_contamination.to_numpy(),obs.raw_pct_counts_mt.to_numpy(),sexmarker])
sensitivity=adjusted.copy()
for stage in ['P0','P7','P14','P21']:
 m=obs.age.eq(stage).to_numpy();c=cov[m];c=(c-c.mean(axis=0))/np.maximum(c.std(axis=0),1e-8);z=adjusted[m];fit=c@np.linalg.lstsq(c,z-z.mean(axis=0),rcond=None)[0];sensitivity[m]=z-fit
representations={'uncorrected':original,'within_stage_harmony':adjusted,'qc_sex_sensitivity':sensitivity}
q=ad.AnnData(sparse.csr_matrix((len(obs),0)),obs=obs.copy());metrics=[];clusterframes=[]
for name,z in representations.items():
 print('GRAPH_UMAP',name,flush=True);q.obsm['X_'+name]=z
 # Independently assess sample mixing within each age; raw cross-age mixing is not a success criterion.
 for stage in ['P0','P7','P14','P21']:
  m=obs.age.eq(stage).to_numpy();s=obs.loc[m,'sample'].astype(str).to_numpy();n=len(s);nn=NearestNeighbors(n_neighbors=31,n_jobs=4).fit(z[m]);neighbors=nn.kneighbors(z[m],return_distance=False)[:,1:]
  cross=(s[:,None]!=s[neighbors]).mean();freq=pd.Series(s).value_counts().to_numpy()/n;expected=(1-(freq**2).sum())*n/(n-1)
  entropy=[]
  for lib in np.unique(s):entropy.append((s[neighbors]==lib).mean(axis=1))
  probs=np.array(entropy);ent=-(probs*np.log(np.maximum(probs,1e-12))).sum(axis=0)/np.log(len(freq))
  metrics.append({'representation':name,'age':stage,'nuclei':n,'different_library_neighbor_fraction':float(cross),'random_mixing_expected_fraction':float(expected),'observed_over_random_mixing':float(cross/expected),'mean_normalized_library_entropy':float(ent.mean())})
 sc.pp.neighbors(q,n_neighbors=15,n_pcs=30,use_rep='X_'+name,random_state=42,key_added=name)
 sc.tl.umap(q,neighbors_key=name,min_dist=.3,maxiter=250,random_state=42);q.obsm['X_umap_'+name]=q.obsm['X_umap'].copy()
 for res in [.4,.6,.8]:
  key=f'{name}_leiden_r{res:g}';sc.tl.leiden(q,neighbors_key=name,resolution=res,flavor='igraph',n_iterations=-1,random_state=42,key_added=key)
  clusterframes.append({'representation':name,'resolution':res,'clusters':int(q.obs[key].nunique())})
q.obs['ventricular_cluster']='V_'+q.obs.within_stage_harmony_leiden_r0_6.astype(str) if 'within_stage_harmony_leiden_r0_6' in q.obs else 'V_'+q.obs['within_stage_harmony_leiden_r0.6'].astype(str)
for name in representations:
 q.obs[name+'_UMAP1']=q.obsm['X_umap_'+name][:,0];q.obs[name+'_UMAP2']=q.obsm['X_umap_'+name][:,1]
pd.DataFrame(metrics).to_csv(OUT/'within_stage_sample_mixing.tsv',sep='\t',index=False);pd.DataFrame(clusterframes).to_csv(OUT/'resolution_sensitivity.tsv',sep='\t',index=False)
ari=[]
for stage in ['P0','P7','P14','P21']:
 m=q.obs.age.eq(stage)
 for other in ['uncorrected','qc_sex_sensitivity']:
  ari.append({'age':stage,'comparison':'within_stage_harmony vs '+other,'adjusted_rand_index':adjusted_rand_score(q.obs.loc[m,'within_stage_harmony_leiden_r0.6'],q.obs.loc[m,other+'_leiden_r0.6'])})
pd.DataFrame(ari).to_csv(OUT/'within_stage_cluster_stability.tsv',sep='\t',index=False)
for col in ['age','sample','subcluster_id','cycling_state']:
 pd.crosstab(q.obs.ventricular_cluster,q.obs[col]).to_csv(OUT/f'ventricular_cluster_by_{col}.tsv',sep='\t')
print('MARKERS',flush=True)
clusters=sorted(q.obs.ventricular_cluster.unique(),key=lambda v:int(v.split('_')[1]));codes=pd.Categorical(q.obs.ventricular_cluster,categories=clusters).codes
membership=sparse.csr_matrix((np.ones(len(q),dtype=np.float64),(codes,np.arange(len(q)))),shape=(len(clusters),len(q)))
ns=np.asarray(membership.sum(axis=1)).ravel();logs=(membership@work.X).toarray();detections=(membership@(x>=1)).toarray()/ns[:,None];sums=(membership@x).toarray();total=logs.sum(axis=0);ranks=[];summaries=[]
for i,c in enumerate(clusters):
 mean=logs[i]/ns[i];rest=(total-logs[i])/(len(q)-ns[i]);diff=mean-rest;eligible=(detections[i]>=.2)&~var.mt.to_numpy()&~var.ribo.to_numpy()&~var.gene_symbol.isin(sexgenes).to_numpy();order=np.flatnonzero(eligible&(diff>.1));order=order[np.argsort(diff[order])[::-1]][:40]
 for rank,j in enumerate(order,1):ranks.append({'ventricular_cluster':c,'rank':rank,'gene_id':var.index[j],'gene_symbol':genes[j],'mean_log1p_10k_difference_vs_rest':float(diff[j]),'fraction_corrected_ge1':float(detections[i,j])})
 sub=q.obs[q.obs.ventricular_cluster.eq(c)];age=sub.age.value_counts().index[0];inage=sub[sub.age.eq(age)];libs=inage['sample'].value_counts();summaries.append({'ventricular_cluster':c,'nuclei':len(sub),'dominant_stage':age,'dominant_stage_fraction':float(len(inage)/len(sub)),'dominant_stage_libraries_ge20':int((libs>=20).sum()),'largest_library_fraction_within_dominant_stage':float(libs.max()/len(inage)),'cycling_marker_positive_fraction':float(sub.cycling_state.eq('cycling-marker positive').mean()),'prior_CAR_counts':';'.join(f'{k}:{v}' for k,v in sub.subcluster_id.value_counts().items() if v>0),'top_descriptive_markers':';'.join(genes[order[:10]])})
pd.DataFrame(ranks).to_csv(OUT/'ventricular_cluster_descriptive_markers.tsv',sep='\t',index=False);pd.DataFrame(summaries).to_csv(OUT/'ventricular_cluster_summary.tsv',sep='\t',index=False)
pb=pd.DataFrame(sums.T,index=var.index,columns=clusters);pb.insert(0,'gene_symbol',genes);pb.index.name='gene_id';pb.to_csv(OUT/'ventricular_cluster_corrected_pseudobulk.tsv.gz',sep='\t',float_format='%.6g')
q.obs.to_csv(OUT/'ventricular_cells_embeddings_and_clusters.csv.gz');q.obsm['X_umap']=q.obsm['X_umap_within_stage_harmony'].copy()
q.uns['analysis_method']={'population':'Previously annotated Ventricular-like only; 43,249 nuclei. Excludes atrial-like and conduction-associated candidates.','primary':'Harmony separately within each stage on sample in a shared ventricular PCA space; restore original stage centroids. Stage is protected and not integrated.','sensitivity':'After primary Harmony, regress measured QC and sex-marker covariates within stage. Can remove biological state; sensitivity only.','counts':'No expression values altered. Expression absent from this metadata-only H5AD; original full count H5AD unchanged.','limitations':'Age and stage-specific batch effects are not identifiable separately. Better sample mixing is not proof of correct biology. New V IDs are computational clusters, not validated subtypes.'}
q.write_h5ad(OUT/'ventricular_embeddings_metadata_only.h5ad',compression='gzip')
print('PLOTS',flush=True)
STAGES={'P0':'#0072B2','P7':'#E69F00','P14':'#009E73','P21':'#CC3311'};palette={c:plt.get_cmap('tab20')(i%20) for i,c in enumerate(clusters)}
def scatter(ax,name,labels,pal,title,ids=False):
 xy=q.obsm['X_umap_'+name];order=np.random.default_rng(42).permutation(len(q));values=np.array(labels.astype(str));ax.scatter(xy[order,0],xy[order,1],s=2,c=[pal[v] for v in values[order]],alpha=.8,linewidths=0,rasterized=True)
 ax.set_title(title,fontsize=14);ax.set_xticks([]);ax.set_yticks([]);ax.set_xlabel('UMAP 1');ax.set_ylabel('UMAP 2');ax.spines[['top','right']].set_visible(False)
 if ids:
  for c in clusters:
   center=np.median(xy[q.obs.ventricular_cluster.eq(c)],axis=0);ax.text(*center,c,fontsize=9,weight='bold',ha='center',bbox=dict(facecolor='white',edgecolor='none',alpha=.9,pad=1))
 return xy
plt.rcParams.update({'font.size':11,'pdf.fonttype':42})
with PdfPages(OUT/'ventricular_stage_batch_report.pdf') as pdf:
 fig,axes=plt.subplots(2,3,figsize=(22,13))
 for col,name in enumerate(representations):
  scatter(axes[0,col],name,q.obs.age,STAGES,name.replace('_',' ')+' — stage')
  sp={s:plt.get_cmap('tab10')(i) for i,s in enumerate(q.obs['sample'].astype(str).unique())};scatter(axes[1,col],name,q.obs['sample'],sp,name.replace('_',' ')+' — library')
 fig.legend(handles=[Line2D([],[],marker='o',linestyle='',color=c,label=a) for a,c in STAGES.items()],loc='upper center',bbox_to_anchor=(.5,.95),ncol=4)
 fig.legend(handles=[Line2D([],[],marker='o',linestyle='',color=c,label=s.replace('GSM','')) for s,c in sp.items()],loc='lower center',bbox_to_anchor=(.5,.015),ncol=4,fontsize=10)
 fig.suptitle('Ventricular-like cardiomyocytes: sample/batch sensitivity',fontsize=22,weight='bold');fig.text(.5,.055,'43,249 nuclei · Harmony adjusts libraries within each stage; ages are not aligned. QC/sex regression is sensitivity only.\nEach column has a separately fitted UMAP; positions and distances cannot be compared between columns.',ha='center',fontsize=11)
 fig.subplots_adjust(top=.89,bottom=.13,hspace=.2,wspace=.17);fig.savefig(OUT/'ventricular_batch_comparison.png',dpi=180);fig.savefig(OUT/'ventricular_batch_comparison.pdf');pdf.savefig(fig);plt.close(fig)
 fig,axes=plt.subplots(1,3,figsize=(23,9));scatter(axes[0],'within_stage_harmony',q.obs.age,STAGES,'Developmental stage');scatter(axes[1],'within_stage_harmony',q.obs.ventricular_cluster,palette,'New ventricular clusters',True);old=sorted(q.obs.subcluster_id.unique(),key=lambda v:int(v.split('_')[1]));oldpal={c:plt.get_cmap('tab10')(i) for i,c in enumerate(old)};xy=scatter(axes[2],'within_stage_harmony',q.obs.subcluster_id,oldpal,'Previous CAR groups on adjusted UMAP')
 for c in old:
  center=np.median(xy[q.obs.subcluster_id.eq(c)],axis=0);axes[2].text(*center,c,fontsize=9,weight='bold',ha='center',bbox=dict(facecolor='white',edgecolor='none',alpha=.9,pad=1))
 fig.legend(handles=[Line2D([],[],marker='o',linestyle='',color=c,label=a) for a,c in STAGES.items()],loc='upper center',bbox_to_anchor=(.5,.94),ncol=4);fig.suptitle('Ventricular-like cardiomyocytes — within-stage sample-adjusted UMAP',fontsize=21,weight='bold');fig.text(.5,.035,'Same adjusted coordinates in all three panels · New V IDs are computational clusters; previous CAR IDs are retained as reference.\nExpression/count matrices unchanged. Stage-specific batch effects remain inseparable from age.',ha='center',fontsize=11);fig.subplots_adjust(top=.84,bottom=.12,wspace=.18);fig.savefig(OUT/'ventricular_adjusted_stage_clusters.png',dpi=180);fig.savefig(OUT/'ventricular_adjusted_stage_clusters.pdf');pdf.savefig(fig);plt.close(fig)
 fig,axes=plt.subplots(2,2,figsize=(17,15),sharex=True,sharey=True);xy=q.obsm['X_umap_within_stage_harmony']
 for ax,stage in zip(axes.flat,STAGES):
  m=q.obs.age.eq(stage).to_numpy();ax.scatter(xy[:,0],xy[:,1],s=1,c='#dddddd',alpha=.3,linewidths=0,rasterized=True);ax.scatter(xy[m,0],xy[m,1],s=3,c=[palette[c] for c in q.obs.loc[m,'ventricular_cluster']],alpha=.85,linewidths=0,rasterized=True)
  for c in clusters:
   cm=m&q.obs.ventricular_cluster.eq(c).to_numpy()
   if cm.sum()>=20:
    center=np.median(xy[cm],axis=0);ax.text(*center,c,fontsize=9,weight='bold',ha='center',bbox=dict(facecolor='white',edgecolor='none',alpha=.9,pad=1))
  ax.set_title(f'{stage} (n={m.sum():,})');ax.set_xticks([]);ax.set_yticks([]);ax.set_xlabel('UMAP 1');ax.set_ylabel('UMAP 2')
 fig.suptitle('Ventricular clusters across stages — same sample-adjusted UMAP',fontsize=19,weight='bold');fig.text(.5,.025,'Color: new ventricular cluster · Gray: other stages · Common coordinates and axes across panels\nNo forced mixing between stages; stage differences include biology and potentially stage-confounded technical effects.',ha='center');fig.tight_layout(rect=(0,.065,1,.95));fig.savefig(OUT/'ventricular_clusters_by_stage.png',dpi=180);fig.savefig(OUT/'ventricular_clusters_by_stage.pdf');pdf.savefig(fig);plt.close(fig)
methods={'input':args.input,'nuclei':len(obs),'embedding_HVG_count':int(hv.sum()),'PCA_components':30,'neighbors':15,'UMAP_min_dist':.3,'UMAP_epochs':250,'seed':42,'Harmony':{'batch':'sample','scope':'separate model per age; common global PCA axes','theta':2,'lambda':1,'nclust':30,'max_iter':20,'postprocessing':'Restore each original age centroid'},'excluded_embedding_genes':sexgenes+['mitochondrial genes','ribosomal genes'],'retained_embedding_biology':'Cell-cycle and developmental genes retained; no age regression.','sensitivity_covariates':['log1p raw UMIs','log1p raw detected genes','DecontX contamination','raw mitochondrial percent','log-normalized Xist,Ddx3y,Kdm5d,Eif2s3y,Uty'],'metric':'30 nearest neighbors within age; observed cross-library neighbor fraction divided by random expectation from library sizes. Higher mixing alone does not establish validity.','packages':{k:importlib.metadata.version(k) for k in ['scanpy','anndata','harmonypy','numpy','scipy','scikit-learn','igraph']}}
(OUT/'methods.json').write_text(json.dumps(methods,indent=2));print('DONE',flush=True)
