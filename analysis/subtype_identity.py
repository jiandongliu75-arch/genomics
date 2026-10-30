"""Evidence-based broad identity review and independent within-lineage subclustering."""
from pathlib import Path
import os,argparse,json,hashlib,tempfile,importlib.metadata
p=argparse.ArgumentParser();p.add_argument('--input',type=Path,required=True);p.add_argument('--metadata',type=Path,required=True);p.add_argument('--criteria',type=Path,required=True);p.add_argument('--outdir',type=Path,required=True);args=p.parse_args()
os.environ.setdefault('MPLCONFIGDIR',str(Path(tempfile.gettempdir())/'genomics-subtype-mpl'))
os.environ.setdefault('NUMBA_CACHE_DIR',str(Path(tempfile.gettempdir())/'genomics-subtype-numba'))
import numpy as np,pandas as pd,h5py,scanpy as sc,anndata as ad
from anndata.io import read_elem,sparse_dataset
from scipy import sparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
cfg=json.loads(args.criteria.read_text());rules=cfg['rules'];out=args.outdir;out.mkdir(parents=True,exist_ok=True);sc.settings.n_jobs=4
meta=sc.read_h5ad(args.metadata);obs=meta.obs.copy()
print('READ_COUNTS',flush=True)
with h5py.File(args.input) as f:
 var=read_elem(f['var']);x=read_elem(f['X']).tocsr().astype(np.float32)
 assert len(obs)==x.shape[0] and obs.index.equals(read_elem(f['obs']).index)
 allgenes=set(sum(cfg['broad_markers'].values(),[])+cfg['cycling_markers']+sum(cfg['cm_marker_states'].values(),[]))
 for defs in cfg['subtypes'].values():
  for d in defs:allgenes.update(d['markers']+d['exclude'])
 ix=np.flatnonzero(var.gene_symbol.isin(allgenes));names=var.iloc[ix].gene_symbol.to_numpy()
 # Canonical marker evidence is also checked in original integer counts.
 raw=sparse_dataset(f['layers/counts_raw'])[:,ix].toarray()
counts=pd.DataFrame(x[:,ix].toarray(),index=obs.index,columns=names).T.groupby(level=0).sum().T
raw=pd.DataFrame(raw,index=obs.index,columns=names).T.groupby(level=0).sum().T
log=np.log1p(counts.to_numpy()*1e4/obs.total_counts.to_numpy()[:,None]);log=pd.DataFrame(log,index=obs.index,columns=counts.columns)
gmean=log.mean();gstd=log.std().clip(lower=.1);det=counts>=rules['marker_detected_count'];rawdet=raw>=1
log.to_csv(out/'curated_marker_log1p_10k.csv.gz',float_format='%.4g',compression={'method':'gzip','compresslevel':1})
missing=sorted(allgenes-set(log.columns));(out/'missing_panel_genes.json').write_text(json.dumps(missing,indent=2))

def panel(d,means):
 gs=[g for g in d['markers'] if g in means.index]
 es=[g for g in d.get('exclude',[]) if g in means.index]
 pos=float(((means[gs]-gmean[gs])/gstd[gs]).mean()) if gs else -999
 neg=float(((means[es]-gmean[es])/gstd[es]).clip(lower=0).mean()) if es else 0
 return pos-.5*neg

def supported(gs,frac,means):return [g for g in gs if g in frac.index and frac[g]>=rules['min_marker_fraction'] and means[g]>=rules['min_marker_mean_log1p_10k']]

# Review each global cluster using independent multi-gene lineage evidence.
broad=[];brows=[];obs['global_cluster']=obs['cluster'].astype(str)
for c,idx in obs.groupby('global_cluster',observed=True).groups.items():
 means=log.loc[idx].mean();frac=det.loc[idx].mean()
 scores={k:panel({'markers':gs},means) for k,gs in cfg['broad_markers'].items()}
 qualifying=[lab for lab in scores if len(supported(cfg['broad_markers'][lab],frac,means))>=2]
 ranking=sorted(qualifying or scores,key=scores.get,reverse=True);label=ranking[0]
 runner=ranking[1] if len(ranking)>1 else 'no other marker-supported lineage'
 margin=scores[label]-scores[runner] if runner in scores else np.nan
 support=supported(cfg['broad_markers'][label],frac,means)
 valid=len(support)>=rules['min_supported_broad_markers'] and (len(ranking)==1 or margin>=rules['min_score_margin']) and scores[label]>=rules['min_positive_score']
 candidate=len(support)>=2 and (len(ranking)==1 or margin>=rules['min_score_margin']) and scores[label]>=rules['min_positive_score']
 assigned=label if valid or candidate else 'Unresolved/mixed'
 status='multi-marker supported' if valid else ('candidate: limited markers' if candidate else 'unresolved')
 broad.append({'global_cluster':str(c),'nuclei':len(idx),'reviewed_identity':assigned,'identity_status':status,'positive_markers':';'.join(support),'top_score':scores[label],'runner_up_identity':runner,'score_margin':margin,'previous_identity':obs.loc[idx,'preliminary_lineage'].mode().iloc[0]})
 for lab,gs in cfg['broad_markers'].items():
  for gene in gs:
   if gene in means.index:brows.append({'global_cluster':str(c),'candidate_lineage':lab,'gene':gene,'mean_log1p_10k':means[gene],'fraction_corrected_ge1':frac[gene],'fraction_raw_ge1':rawdet.loc[idx,gene].mean()})
broad=pd.DataFrame(broad);broad.to_csv(out/'reviewed_broad_identities.tsv',sep='\t',index=False);pd.DataFrame(brows).to_csv(out/'broad_marker_evidence.tsv',sep='\t',index=False)
obs['reviewed_identity']=obs.global_cluster.map(broad.set_index('global_cluster').reviewed_identity)
obs['broad_identity_status']=obs.global_cluster.map(broad.set_index('global_cluster').identity_status)
obs['subcluster_parent']=obs.reviewed_identity.replace({'Pericyte':'Mural','Smooth muscle':'Mural'})
obs['cycling_marker_count']=det[[g for g in cfg['cycling_markers'] if g in det]].sum(axis=1)
obs['cycling_state']=np.where(obs.cycling_marker_count>=2,'cycling-marker positive','not classified cycling')
obs['subcluster_id']='not_subclustered';obs['proposed_subtype']=obs.reviewed_identity+' (not subclustered)';obs['subtype_confidence']='broad identity only';obs['cm_marker_state']='not applicable'
obs['subtype_review_flag']=True
summaries=[];evidence=[];markers_out=[];panelrows=[]
parents=['Cardiomyocyte','Fibroblast','Endothelial','Mural','Immune']
for parent in parents:
 mask=obs.subcluster_parent.eq(parent).to_numpy();indices=np.flatnonzero(mask);n=len(indices)
 print('SUBCLUSTER',parent,n,flush=True)
 if n<100:continue
 path=out/parent.replace(' ','_');path.mkdir(exist_ok=True)
 work=ad.AnnData(x[indices].copy(),obs=obs.iloc[indices][['sample','age']].copy(),var=var[['gene_symbol','mt','ribo']].copy())
 sc.pp.normalize_total(work,target_sum=1e4);sc.pp.log1p(work)
 sc.pp.highly_variable_genes(work,n_top_genes=min(2000,work.n_vars),flavor='seurat',batch_key='sample')
 hv=work.var.highly_variable.to_numpy() & ~work.var.mt.to_numpy()
 tmp=work[:,hv].copy();sc.pp.scale(tmp,zero_center=False,max_value=10)
 sc.pp.pca(tmp,n_comps=min(30,n-2),random_state=42,svd_solver='arpack')
 sc.pp.neighbors(tmp,n_neighbors=min(15,n-1),n_pcs=min(30,n-2),random_state=42)
 for res in [.3,.6,1.0]:sc.tl.leiden(tmp,resolution=res,flavor='igraph',n_iterations=-1,random_state=42,key_added=f'leiden_r{res:g}')
 sc.tl.umap(tmp,min_dist=.3,maxiter=200,random_state=42)
 primary=tmp.obs['leiden_r0.6'];cluster_names=sorted(primary.unique(),key=int)
 ids=parent[:3].upper()+'_'+primary.astype(str)
 obs.loc[work.obs_names,'subcluster_id']=ids.to_numpy()
 obs.loc[work.obs_names,'subcluster_UMAP1']=tmp.obsm['X_umap'][:,0];obs.loc[work.obs_names,'subcluster_UMAP2']=tmp.obsm['X_umap'][:,1]
 group=sparse.csr_matrix((np.ones(n,dtype=np.float64),(pd.Categorical(primary,categories=cluster_names).codes,np.arange(n))),shape=(len(cluster_names),n))
 ns=np.asarray(group.sum(axis=1)).ravel();logsum=(group@work.X).toarray();totallog=np.asarray(work.X.astype(np.float64).sum(axis=0)).ravel()
 detected=x[indices].copy();detected.data=(detected.data>=1).astype(np.float32);detected.eliminate_zeros();fractions=(group@detected).toarray()/ns[:,None];totaldet=np.asarray(detected.sum(axis=0)).ravel();del detected
 sums=(group@x[indices]).toarray();table=pd.DataFrame(sums.T,index=var.index,columns=[parent[:3].upper()+'_'+str(c) for c in cluster_names]);table.insert(0,'gene_symbol',var.gene_symbol.to_numpy());table.to_csv(path/'subcluster_pseudobulk_corrected_counts.tsv.gz',sep='\t',float_format='%.6g',compression={'method':'gzip','compresslevel':1})
 for i,c in enumerate(cluster_names):
  cid=parent[:3].upper()+'_'+str(c);cellids=work.obs_names[primary.eq(c).to_numpy()];others=work.obs_names[~primary.eq(c).to_numpy()]
  means=log.loc[cellids].mean();frac=det.loc[cellids].mean();rawfrac=rawdet.loc[cellids].mean();samplecounts=obs.loc[cellids,'sample'].value_counts();nsamples=int((samplecounts>=10).sum());maxsample=float(samplecounts.max()/len(cellids));enrich=means-log.loc[others].mean()
  candidates=[]
  for d in cfg['subtypes'][parent]:
   supports=supported(d['markers'],frac,means);anchors=supported(d['anchors'],frac,means);gs=[g for g in d['markers'] if g in log]
   contrasts=[]
   for sample,number in samplecounts.items():
    if number<10:continue
    same=obs.loc[others,'sample'].eq(sample);control=others[same.to_numpy()]
    if len(control)<10:continue
    value=float(log.loc[cellids[obs.loc[cellids,'sample'].eq(sample).to_numpy()],gs].mean().mean()-log.loc[control,gs].mean().mean());contrasts.append((sample,value))
   positive=sum(value>0 for _,value in contrasts)
   enrichvalue=float(enrich[gs].mean())
   marker_pass=len(supports)>=rules['min_supported_panel_markers'] and len(anchors)>=d['anchor_min']
   score=panel(d,means)
   candidates.append({'definition':d,'score':score,'supported':supports,'anchors':anchors,'positive_contrasts':positive,'enrichment':enrichvalue,'marker_pass':marker_pass})
   panelrows.append({'parent':parent,'subcluster_id':cid,'candidate_subtype':d['name'],'panel_score':score,'supported_markers':';'.join(supports),'supported_anchors':';'.join(anchors),'marker_rule_pass':marker_pass,'mean_panel_log_difference':enrichvalue,'libraries_with_positive_within_library_contrast':positive,'library_contrasts':';'.join(f'{s}:{v:.4g}' for s,v in contrasts)})
   for gene in set(d['markers']+d['exclude']):
    if gene in means.index:evidence.append({'parent':parent,'subcluster_id':cid,'candidate_subtype':d['name'],'gene':gene,'role':'positive' if gene in d['markers'] else 'exclusion/check','mean_log1p_10k':means[gene],'fraction_corrected_ge1':frac[gene],'fraction_raw_ge1':rawfrac[gene],'mean_log_difference_vs_parent_rest':enrich[gene]})
  candidates.sort(key=lambda z:z['score'],reverse=True);eligible=[z for z in candidates if z['marker_pass']]
  chosen=eligible[0] if eligible else None
  margin=(chosen['score']-max([v['score'] for v in candidates if v is not chosen],default=-999)) if chosen else np.nan
  if chosen:
   label=chosen['definition']['name'];quality=len(cellids)>=rules['min_cells_for_supported_subtype'] and nsamples>=rules['min_samples_with_10_cells'] and maxsample<=rules['max_single_sample_fraction']
   distinct=margin>=rules['min_score_margin'] and chosen['score']>=rules['min_positive_score'] and chosen['enrichment']>=rules['min_relative_panel_enrichment'] and chosen['positive_contrasts']>=rules['min_within_sample_positive_contrasts']
   confidence='provisional-supported' if quality and distinct else 'candidate/review'
  else:label=parent+' unresolved';confidence='unresolved'
  state=[]
  if parent=='Cardiomyocyte':
   for name,gs in cfg['cm_marker_states'].items():
    sup=supported(gs,frac,means);anchor=gs[0]
    if anchor in sup and len(sup)>=2 and enrich.get(anchor,0)>.1:state.append(name)
  cycling=float(obs.loc[cellids,'cycling_state'].eq('cycling-marker positive').mean())
  obs.loc[cellids,'proposed_subtype']=label;obs.loc[cellids,'subtype_confidence']=confidence;obs.loc[cellids,'subtype_review_flag']=confidence!='provisional-supported' or 'candidate' in label.lower();obs.loc[cellids,'cm_marker_state']=';'.join(state) if state else 'no assigned marker state'
  summary={'parent':parent,'subcluster_id':cid,'nuclei':len(cellids),'proposed_subtype':label,'confidence':confidence,'panel_score':chosen['score'] if chosen else np.nan,'score_margin':margin,'supported_markers':';'.join(chosen['supported']) if chosen else '', 'supported_anchors':';'.join(chosen['anchors']) if chosen else '', 'positive_within_library_contrasts':chosen['positive_contrasts'] if chosen else 0,'panel_enrichment':chosen['enrichment'] if chosen else np.nan,'libraries_with_ge10_nuclei':nsamples,'largest_library_fraction':maxsample,'cycling_positive_fraction':cycling,'cycling_enriched':cycling>=.1,'cm_marker_state':';'.join(state),'dominant_age':obs.loc[cellids,'age'].value_counts().index[0],'dominant_age_fraction':obs.loc[cellids,'age'].value_counts().max()/len(cellids),'median_ambient_estimate':obs.loc[cellids,'decontx_contamination'].median()}
  summaries.append(summary)
  mean=logsum[i]/ns[i];rest=(totallog-logsum[i])/(n-ns[i]);diff=mean-rest;frest=(totaldet-fractions[i]*ns[i])/(n-ns[i]);eligiblegenes=(fractions[i]>=.2)&((fractions[i]-frest)>=.1)&(~var.mt.to_numpy())&(~var.ribo.to_numpy());order=np.flatnonzero(eligiblegenes);order=order[np.argsort(diff[order])[::-1]][:40]
  for rank,j in enumerate(order,1):markers_out.append({'parent':parent,'subcluster_id':cid,'rank':rank,'gene_id':var.index[j],'gene_symbol':var.gene_symbol.iloc[j],'mean_log1p_10k_difference':diff[j],'fraction_ge1_subcluster':fractions[i,j],'fraction_ge1_parent_rest':frest[j]})
 # Save per-parent metadata and independently computed PCA/UMAP/graph.
 submeta=ad.AnnData(sparse.csr_matrix((n,0)),obs=obs.loc[work.obs_names].copy());submeta.obsm['X_umap']=tmp.obsm['X_umap'];submeta.obsm['X_pca']=tmp.obsm['X_pca'];submeta.obsp.update(tmp.obsp);submeta.uns['neighbors']=dict(tmp.uns['neighbors'])
 for res in [.3,.6,1.0]:submeta.obs[f'leiden_r{res:g}']=tmp.obs[f'leiden_r{res:g}'].to_numpy()
 submeta.uns['method']={'primary_resolution':.6,'seed':42,'integration':'none','parent':parent,'expression_matrix':'metadata only; original and corrected expression not stored here'}
 submeta.write_h5ad(path/'subclusters_metadata_only.h5ad',compression='gzip');submeta.obs.to_csv(path/'cell_subcluster_annotations.csv.gz',compression={'method':'gzip','compresslevel':1})
 for column in ['age','sample','proposed_subtype','cycling_state']:
  pd.crosstab(submeta.obs.subcluster_id,submeta.obs[column]).to_csv(path/f'subcluster_by_{column}.tsv',sep='\t')
 fig,axes=plt.subplots(1,3,figsize=(20,6));xy=tmp.obsm['X_umap'];cmap=plt.get_cmap('tab20')
 for j,c in enumerate(cluster_names):
  cm=primary.eq(c).to_numpy();axes[0].scatter(xy[cm,0],xy[cm,1],s=2,color=cmap(j%20),rasterized=True);center=np.median(xy[cm],axis=0);axes[0].text(*center,str(c),ha='center',fontsize=8,bbox={'facecolor':'white','alpha':.8,'pad':1,'edgecolor':'none'})
 labels=list(submeta.obs.proposed_subtype.unique())
 for j,label in enumerate(labels):
  lm=submeta.obs.proposed_subtype.eq(label).to_numpy();axes[1].scatter(xy[lm,0],xy[lm,1],s=2,label=label,color=cmap(j%20),rasterized=True)
 axes[1].legend(fontsize=6,markerscale=4,loc='best')
 for age in ['P0','P7','P14','P21']:
  am=submeta.obs.age.eq(age).to_numpy();axes[2].scatter(xy[am,0],xy[am,1],s=2,label=age,rasterized=True)
 axes[2].legend(markerscale=4,fontsize=7)
 for ax,title in zip(axes,['Within-lineage subclusters','Proposed subtype (including candidates)','Age, no integration']):ax.set(title=title,xlabel='UMAP1',ylabel='UMAP2')
 fig.suptitle(parent);fig.tight_layout();fig.savefig(path/'subcluster_umap.png',dpi=160);fig.savefig(path/'subcluster_umap.pdf');plt.close(fig)
 # Diagnostic dot plot: size = fraction >=1 corrected count; color = within-gene scaled mean.
 plotgenes=list(dict.fromkeys(g for d in cfg['subtypes'][parent] for g in d['markers'] if g in counts))
 if parent=='Cardiomyocyte':plotgenes+=list(dict.fromkeys(g for gs in cfg['cm_marker_states'].values() for g in gs if g in counts and g not in plotgenes))
 values=np.vstack([log.loc[work.obs_names[primary.eq(c).to_numpy()],plotgenes].mean().to_numpy() for c in cluster_names]);proportions=np.vstack([det.loc[work.obs_names[primary.eq(c).to_numpy()],plotgenes].mean().to_numpy() for c in cluster_names]);scaled=values/np.maximum(values.max(axis=0),.01)
 fig,ax=plt.subplots(figsize=(max(12,len(plotgenes)*.30),max(4,len(cluster_names)*.32)))
 for i in range(len(cluster_names)):
  im=ax.scatter(np.arange(len(plotgenes)),np.repeat(i,len(plotgenes)),s=proportions[i]*90,c=scaled[i],cmap='viridis',vmin=0,vmax=1)
 ax.set_xticks(range(len(plotgenes)),plotgenes,rotation=90,fontsize=7);ax.set_yticks(range(len(cluster_names)),[parent[:3].upper()+'_'+str(c) for c in cluster_names],fontsize=8);ax.invert_yaxis();ax.set_title(parent+' | size: fraction >=1 count; color: within-gene scaled mean');fig.colorbar(im,ax=ax,shrink=.5);fig.tight_layout();fig.savefig(path/'subtype_marker_dotplot.png',dpi=160);fig.savefig(path/'subtype_marker_dotplot.pdf');plt.close(fig)
 print('DONE_PARENT',parent,len(cluster_names),flush=True)
 del work,tmp,submeta,logsum,fractions,sums
# Main annotations include every retained nucleus, including minor/unresolved populations.
obs.to_csv(out/'all_cell_identity_and_subtype_annotations.csv.gz',compression={'method':'gzip','compresslevel':1})
summary=pd.DataFrame(summaries);summary.to_csv(out/'subcluster_identity_evidence.tsv',sep='\t',index=False);pd.DataFrame(evidence).to_csv(out/'subtype_marker_evidence.tsv',sep='\t',index=False);pd.DataFrame(panelrows).to_csv(out/'candidate_panel_comparisons.tsv',sep='\t',index=False);pd.DataFrame(markers_out).to_csv(out/'subcluster_descriptive_markers.tsv',sep='\t',index=False)
for col in ['reviewed_identity','proposed_subtype','subtype_confidence']:
 pd.crosstab(obs[col],obs['sample']).to_csv(out/f'{col}_counts_by_sample.tsv',sep='\t')
allmeta=ad.AnnData(sparse.csr_matrix((len(obs),0)),obs=obs.copy());allmeta.obsm['X_umap']=meta.obsm['X_umap'];allmeta.uns['annotation_method']={'criteria_file':'subtype_criteria.json','labels':'provisional: inspect subtype_confidence and subtype_review_flag','expression_matrix':'absent in metadata-only H5AD'};allmeta.write_h5ad(out/'all_cell_annotations_metadata_only.h5ad',compression='gzip')
fig,axes=plt.subplots(1,2,figsize=(16,7));xy=meta.obsm['X_umap'];cmap=plt.get_cmap('tab20')
for i,label in enumerate(obs.reviewed_identity.unique()):
 mask=obs.reviewed_identity.eq(label).to_numpy();axes[0].scatter(xy[mask,0],xy[mask,1],s=1,label=label,color=cmap(i%20),rasterized=True)
axes[0].legend(fontsize=7,markerscale=5);axes[0].set_title('Reviewed broad identities')
for label,col in [('provisional-supported','#27867c'),('candidate/review','#e7a34a'),('unresolved','#b95e60'),('broad identity only','#bdbdbd')]:
 mask=obs.subtype_confidence.eq(label).to_numpy();axes[1].scatter(xy[mask,0],xy[mask,1],s=1,label=label,color=col,rasterized=True)
axes[1].legend(fontsize=7,markerscale=5);axes[1].set_title('Subtype evidence status, not a probability')
for ax in axes:ax.set(xlabel='Original global UMAP1',ylabel='Original global UMAP2')
fig.tight_layout();fig.savefig(out/'reviewed_identities_and_confidence.png',dpi=160);fig.savefig(out/'reviewed_identities_and_confidence.pdf');plt.close(fig)
# Validate that no nucleus was lost or duplicated, and all subcluster counts reconcile.
assert obs.index.is_unique and len(obs)==73052 and obs.reviewed_identity.notna().all()
for parent in parents:
 cells=obs.subcluster_parent.eq(parent);rows=summary.parent.eq(parent)
 assert int(cells.sum())==int(summary.loc[rows,'nuclei'].sum())
assert set(obs.index)==set(meta.obs_names)
roundtrip=pd.read_csv(out/'all_cell_identity_and_subtype_annotations.csv.gz',index_col=0);assert roundtrip.index.equals(obs.index)
valid=summary.confidence.eq('provisional-supported')
assert (summary.loc[valid,'nuclei']>=50).all() and (summary.loc[valid,'libraries_with_ge10_nuclei']>=2).all() and (summary.loc[valid,'positive_within_library_contrasts']>=2).all()
versions={name:importlib.metadata.version(name) for name in ['scanpy','anndata','numpy','pandas','scipy','matplotlib','igraph','leidenalg']}
(out/'requirements.txt').write_text('\n'.join(k+'=='+v for k,v in versions.items())+'\n')
def sha(path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  for chunk in iter(lambda:f.read(1048576),b''):h.update(chunk)
 return h.hexdigest()
manifest={'status':'validated','nuclei':len(obs),'subclusters':len(summary),'confidence_counts':obs.subtype_confidence.value_counts().to_dict(),'broad_identity_counts':obs.reviewed_identity.value_counts().to_dict(),'parameters':{'primary_resolution':.6,'comparison_resolutions':[.3,.6,1.0],'hvg':2000,'pca':30,'neighbors':15,'seed':42,'integration':'none'},'input_sha256':sha(args.input),'criteria_sha256':sha(args.criteria),'versions':versions,'files':[]}
for path in sorted(out.rglob('*')):
 if path.is_file() and path.name!='manifest.json':manifest['files'].append({'file':str(path.relative_to(out)),'bytes':path.stat().st_size,'sha256':sha(path)})
assert all(f['bytes']<90_000_000 for f in manifest['files'])
(out/'manifest.json').write_text(json.dumps(manifest,indent=2))
print('COMPLETE',json.dumps({k:v for k,v in manifest.items() if k in ['nuclei','subclusters','confidence_counts','broad_identity_counts']}),flush=True)
print(summary[['parent','subcluster_id','nuclei','proposed_subtype','confidence','supported_markers','positive_within_library_contrasts']].to_string(index=False),flush=True)
