"""Independent within-stage ventricular clustering with replicate and nuisance sensitivity."""
from pathlib import Path
import os,argparse,json,importlib.metadata,re
os.environ.setdefault('MPLCONFIGDIR','/tmp/heart-mpl');os.environ.setdefault('NUMBA_CACHE_DIR','/tmp/heart-numba')
import numpy as np,pandas as pd,h5py,scanpy as sc,anndata as ad,harmonypy
from anndata.io import read_elem,sparse_dataset
from scipy import sparse
from sklearn.neighbors import NearestNeighbors
from sklearn.metrics import adjusted_rand_score
p=argparse.ArgumentParser();p.add_argument('--input',default='/workspace/heart-atlas/results/cell_identities_subtypes/snrna_qc_annotated_subclustered.h5ad');args=p.parse_args()
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/ventricular_within_stage';OUT.mkdir(exist_ok=True);sc.settings.n_jobs=4
SEX=['Xist','Tsix','Ddx3y','Kdm5d','Eif2s3y','Uty','Usp9y','Zfy1','Zfy2','Sry','Rbmy']
PANELS={'Cycling/DNA-replication program':['Mki67','Top2a','Cdk1','Ccnb1','Cenpf','Mcm2','Mcm5','Pcna'], 'Lipid/substrate-handling program':['Cd36','Slc27a1','Cpt1b','Pdk4','Hmgcs2','Acadm','Hadha'], 'Calcium/contractile program':['Ryr2','Atp2a2','Casq2','Pln','Jph2','Srl','Actn2','Tnni3'], 'Adhesion/structural program':['Ncam1','Vcl','Smpx','Bves','Tecrl','Fhl2'], 'Neonatal contractile program':['Tnni1','Myh7','Nppa','Bmp10'], 'Response-associated program':['Ankrd1','Nppb','Fos','Jun','Junb','Atf3','Egr1','Hspa1a'], 'Compact-associated candidate':['Hey2','Mycn'], 'Trabecular-associated candidate':['Nppa','Bmp10','Erbb4'], 'Conduction-overlap check':['Etv1','Cntn2','Gja5','Scn5a']}
with h5py.File(args.input) as f:
 obs=read_elem(f['obs']);var=read_elem(f['var']);ix=np.flatnonzero(obs.proposed_subtype.eq('Ventricular-like').to_numpy());obs=obs.iloc[ix].copy();print('READ',len(obs),flush=True);x=sparse_dataset(f['X'])[ix,:].tocsr().astype(np.float32)
assert len(obs)==43249
w=ad.AnnData(x.copy(),obs=obs[['sample','age']].copy(),var=var.copy());sc.pp.normalize_total(w,target_sum=1e4);sc.pp.log1p(w)
genes=var.gene_symbol.astype(str).to_numpy();lookup={g:np.flatnonzero(genes==g) for g in set(sum(PANELS.values(),[])+SEX+['Col1a1','Pdgfra','Pecam1','Ptprc','Tnnt2','Myl2','Myl3'])}
allann=[];summaries=[];markers=[];profiles=[];mix=[];resolution=[];panelrows=[];hvgrows=[];globalstability=[]
def mean(z,m):return np.asarray(z[m].mean(axis=0)).ravel()
def best_jaccard(ref,other,label):
 a=ref==label;best=0;match=''
 for b in np.unique(other):
  bm=other==b;j=np.sum(a&bm)/np.sum(a|bm)
  if j>best:best=float(j);match=str(b)
 return best,match
for stage in ['P0','P7','P14','P21']:
 print('STAGE',stage,flush=True);ind=np.flatnonzero(obs.age.eq(stage).to_numpy());o=obs.iloc[ind].copy()
 for col in o.select_dtypes(include='category').columns:o[col]=o[col].cat.remove_unused_categories()
 stageout=OUT/stage;stageout.mkdir(exist_ok=True)
 zlog=w.X[ind].tocsr();counts=x[ind].tocsr();swork=ad.AnnData(zlog.copy(),obs=o[['sample','age']].copy(),var=var.copy())
 sc.pp.highly_variable_genes(swork,n_top_genes=2000,flavor='seurat',batch_key='sample');hv=swork.var.highly_variable.to_numpy()&~var.mt.to_numpy()&~var.ribo.to_numpy()&~var.gene_symbol.isin(SEX).to_numpy()
 for j in np.flatnonzero(hv):hvgrows.append({'stage':stage,'gene_id':var.index[j],'gene_symbol':genes[j]})
 tmp=swork[:,hv].copy();sc.pp.scale(tmp,zero_center=False,max_value=10);sc.pp.pca(tmp,n_comps=30,random_state=42,svd_solver='arpack');rawpc=tmp.obsm['X_pca'].copy()
 print('HARMONY',stage,int(hv.sum()),flush=True);h=harmonypy.run_harmony(rawpc,o[['sample']],['sample'],theta=2,lamb=1,nclust=30,max_iter_harmony=20,max_iter_kmeans=10,random_state=42,ncores=4,verbose=False);pc=np.asarray(h.Z_corr);pc=pc if pc.shape==rawpc.shape else pc.T;pc=pc-pc.mean(axis=0)+rawpc.mean(axis=0)
 cov=[np.log1p(o.raw_total_counts.to_numpy()),np.log1p(o.raw_n_genes_by_counts.to_numpy()),o.decontx_contamination.to_numpy(),o.raw_pct_counts_mt.to_numpy()]
 for gene in ['Xist','Ddx3y','Kdm5d','Eif2s3y','Uty']:cov.append(np.asarray(zlog[:,lookup[gene]].sum(axis=1)).ravel())
 cov=np.array(cov).T;cov=(cov-cov.mean(axis=0))/np.maximum(cov.std(axis=0),1e-8);sensitive=pc-cov@np.linalg.lstsq(cov,pc-pc.mean(axis=0),rcond=None)[0]
 q=ad.AnnData(sparse.csr_matrix((len(o),0)),obs=o.copy());rep={'uncorrected':rawpc,'harmony':pc,'qc_sex_sensitivity':sensitive}
 print('GRAPHS',stage,flush=True)
 for name,v in rep.items():
  q.obsm['X_'+name]=v;sc.pp.neighbors(q,n_neighbors=15,n_pcs=30,use_rep='X_'+name,random_state=42,key_added=name)
  for res in [.3,.6,.9,1.2]:
   key=f'{name}_r{res:g}';sc.tl.leiden(q,neighbors_key=name,resolution=res,flavor='igraph',n_iterations=-1,random_state=42,key_added=key);resolution.append({'stage':stage,'representation':name,'resolution':res,'clusters':int(q.obs[key].nunique())})
  sample=o['sample'].astype(str).to_numpy();nn=NearestNeighbors(n_neighbors=31,n_jobs=4).fit(v).kneighbors(v,return_distance=False)[:,1:];cross=np.mean(sample[:,None]!=sample[nn]);freq=o['sample'].value_counts().to_numpy()/len(o);expected=(1-(freq**2).sum())*len(o)/(len(o)-1);mix.append({'stage':stage,'representation':name,'different_library_neighbor_fraction':float(cross),'random_mixing_expected':float(expected),'observed_over_random':float(cross/expected)})
 for seed in [7,123]:sc.tl.leiden(q,neighbors_key='harmony',resolution=.6,flavor='igraph',n_iterations=-1,random_state=seed,key_added=f'harmony_seed{seed}')
 sc.tl.umap(q,neighbors_key='harmony',min_dist=.3,maxiter=250,random_state=42);q.obs['stage_subcluster_id']=stage+'_V'+q.obs['harmony_r0.6'].astype(str)
 ref=q.obs['harmony_r0.6'].astype(str).to_numpy();alllabels={key:q.obs[key].astype(str).to_numpy() for key in ['uncorrected_r0.6','qc_sex_sensitivity_r0.6','harmony_seed7','harmony_seed123','harmony_r0.3','harmony_r0.9','harmony_r1.2']}
 for key,v in alllabels.items():globalstability.append({'stage':stage,'comparison':'primary_vs_'+key,'adjusted_rand_index':adjusted_rand_score(ref,v)})
 q.obs['within_stage_UMAP1']=q.obsm['X_umap'][:,0];q.obs['within_stage_UMAP2']=q.obsm['X_umap'][:,1]
 clusterids=sorted(q.obs.stage_subcluster_id.unique(),key=lambda k:int(k.split('V')[-1]));codes=pd.Categorical(q.obs.stage_subcluster_id,categories=clusterids).codes;g=sparse.csr_matrix((np.ones(len(o),dtype=np.float64),(codes,np.arange(len(o)))),shape=(len(clusterids),len(o)));ns=np.asarray(g.sum(axis=1)).ravel();logs=(g@zlog).toarray();sums=(g@counts).toarray();fractions=(g@(counts>=1)).toarray()/ns[:,None];totallog=logs.sum(axis=0)
 means=logs/ns[:,None];std=np.maximum(np.std(means,axis=0),.1)
 libdata={}
 for lib in o['sample'].unique():
  lm=o['sample'].eq(lib).to_numpy();lns=np.asarray(g[:,lm].sum(axis=1)).ravel();llog=(g[:,lm]@zlog[lm]).toarray();libdata[lib]=(lns,llog)
 for k,cid in enumerate(clusterids):
  m=q.obs.stage_subcluster_id.eq(cid).to_numpy();sub=q.obs.loc[m];diff=means[k]-(totallog-logs[k])/(len(o)-ns[k]);restfrac=(fractions*ns[:,None]).sum(axis=0)-fractions[k]*ns[k];restfrac/=len(o)-ns[k]
  contrasts={}
  for lib,(lns,llog) in libdata.items():
   na=lns[k];nb=lns.sum()-na
   if min(na,nb)>=30:contrasts[lib]=llog[k]/na-(llog.sum(axis=0)-llog[k])/nb
  support=np.sum(np.array(list(contrasts.values()))>.05,axis=0) if contrasts else np.zeros(len(var),dtype=int)
  eligible=(diff>.15)&(fractions[k]>=.2)&((fractions[k]-restfrac)>=.08)&~var.mt.to_numpy()&~var.ribo.to_numpy()&~var.gene_symbol.isin(SEX).to_numpy()
  order=np.flatnonzero(eligible);order=order[np.argsort(diff[order])[::-1]][:60]
  coherent=[j for j in order if support[j]>=2]
  for rank,j in enumerate(order,1):
   row={'stage':stage,'stage_subcluster_id':cid,'rank':rank,'gene_id':var.index[j],'gene_symbol':genes[j],'mean_log1p_10k_difference_vs_stage_rest':diff[j],'fraction_corrected_ge1_cluster':fractions[k,j],'fraction_corrected_ge1_stage_rest':restfrac[j],'eligible_libraries':len(contrasts),'positive_libraries':int(support[j])}
   for lib,d in contrasts.items():row[lib+'_difference']=d[j]
   markers.append(row)
  panelinfo=[]
  for name,panel in PANELS.items():
   js=[int(j) for gene in panel for j in lookup[gene]];positive=[j for j in js if diff[j]>.1 and fractions[k,j]>=.1 and support[j]>=2];score=float(np.mean(diff[js]/std[js])) if js else -np.inf
   if name=='Compact-associated candidate':qualified=len(positive)==2
   elif name=='Trabecular-associated candidate':qualified=len(positive)>=2 and any(genes[j] in ['Bmp10','Erbb4'] for j in positive)
   else:qualified=len(positive)>=2
   qualified=qualified and score>0
   panelinfo.append((name,score,qualified,positive));panelrows.append({'stage':stage,'stage_subcluster_id':cid,'program':name,'relative_score':score,'supported_panel_genes':';'.join(genes[positive]),'panel_rule_pass':qualified})
  cycling=float(sub.cycling_state.eq('cycling-marker positive').mean())
  qualified=sorted([v for v in panelinfo if v[2] and v[0] not in ['Cycling/DNA-replication program','Conduction-overlap check']],key=lambda v:v[1],reverse=True)
  if cycling>=.1 and next(v for v in panelinfo if v[0]=='Cycling/DNA-replication program')[2]:label='Cycling/DNA-replication program'
  elif qualified:label=qualified[0][0]
  else:
   labelgenes=[genes[j] for j in coherent if not re.match(r'^(Gm\d|\d)|Rik$',genes[j])] or [genes[j] for j in order if not re.match(r'^(Gm\d|\d)|Rik$',genes[j])]
   label='Marker-defined: '+('/'.join(labelgenes[:2])) if labelgenes else 'Uncharacterized ventricular state'
  selected_program=label if label in PANELS else ''
  selected_genes=set(genes[next(v for v in panelinfo if v[0]==label)[3]]) if selected_program else set()
  if label=='Calcium/contractile program' and selected_genes=={'Srl','Actn2'}:label='Srl/Actn2-enriched structural program'
  elif label=='Adhesion/structural program' and len(selected_genes)<=2 and {'Fgf14','Thsd4'}.issubset(set(genes[coherent])):label='Fgf14/Thsd4-enriched state'
  elif label=='Neonatal contractile program' and {'Myh7','Nppa'}.issubset(selected_genes):label='Myh7/Nppa-enriched program'
  matches={key:best_jaccard(ref,v,str(q.obs.loc[m,'harmony_r0.6'].iloc[0])) for key,v in alllabels.items()};seedvals=[matches['harmony_seed7'][0],matches['harmony_seed123'][0]];libs=sub['sample'].value_counts();enough=len(sub)>=100 and (libs>=30).sum()==2 and libs.max()/len(sub)<=.9;stable=min(seedvals)>=.5 and np.mean(seedvals)>=.65;markerpass=len(coherent)>=3
  replicate=enough and stable and markerpass;nuisance=matches['qc_sex_sensitivity_r0.6'][0]<.5;rawsensitive=matches['uncorrected_r0.6'][0]<.4
  status='supported transcriptional state' if replicate and not nuisance and not rawsensitive else 'replicate-supported; adjustment-sensitive' if replicate else 'candidate/review'
  if 'candidate' in label.lower():status='candidate anatomical association; '+status
  reasons=[]
  if not enough:reasons.append('size/library representation')
  if not stable:reasons.append('seed sensitivity')
  if not markerpass:reasons.append('fewer than 3 cross-library markers')
  if nuisance:reasons.append('QC/sex-sensitive')
  if rawsensitive:reasons.append('sample-adjustment-sensitive')
  row={'stage':stage,'stage_subcluster_id':cid,'nuclei':len(sub),'stage_state_label':label,'selected_program_panel':selected_program,'stage_support_status':status,'replicate_program_support':bool(replicate),'qc_sex_sensitive':bool(nuisance),'sample_adjustment_sensitive':bool(rawsensitive),'review_reasons':';'.join(reasons),'libraries_ge30':int((libs>=30).sum()),'largest_library_fraction':float(libs.max()/len(sub)),'both_library_marker_count':len(coherent),'both_library_markers':';'.join(genes[coherent[:15]]),'cycling_positive_fraction':cycling,'seed_mean_best_jaccard':float(np.mean(seedvals)),'seed_min_best_jaccard':float(min(seedvals)),'median_raw_counts':float(sub.raw_total_counts.median()),'median_raw_genes':float(sub.raw_n_genes_by_counts.median()),'median_contamination':float(sub.decontx_contamination.median())}
  for key,(val,match) in matches.items():row[key+'_best_jaccard']=val;row[key+'_best_match']=match
  summaries.append(row);q.obs.loc[m,'stage_state_label']=label;q.obs.loc[m,'stage_support_status']=status
  for gene,js in lookup.items():
   for j in js:profiles.append({'stage':stage,'stage_subcluster_id':cid,'gene':gene,'mean_log1p_10k':means[k,j],'fraction_corrected_ge1':fractions[k,j],'mean_log_difference_vs_stage_rest':diff[j],'positive_libraries':int(support[j])})
 pb=pd.DataFrame(sums.T,index=var.index,columns=clusterids);pb.insert(0,'gene_symbol',genes);pb.index.name='gene_id';pb.to_csv(stageout/'corrected_pseudobulk.tsv.gz',sep='\t',float_format='%.6g')
 pd.crosstab(q.obs.stage_subcluster_id,q.obs['sample']).to_csv(stageout/'cluster_by_sample.tsv',sep='\t');pd.crosstab(q.obs.stage_subcluster_id,q.obs.subcluster_id).to_csv(stageout/'cluster_by_prior_CAR.tsv',sep='\t')
 q.obs.to_csv(stageout/'cell_annotations.csv.gz');q.uns['within_stage_method']={'stage':stage,'genes':int(hv.sum()),'primary':'Independent stage HVG/PCA; Harmony by sample; Leiden r0.6 seed42; stage UMAP','expression':'absent in metadata-only object; expression/count values unchanged'};q.write_h5ad(stageout/'embeddings_metadata_only.h5ad',compression='gzip');allann.append(q.obs);print('FINISHED',stage,len(clusterids),flush=True)
pd.concat(allann).to_csv(OUT/'cell_annotations.csv.gz');pd.DataFrame(summaries).to_csv(OUT/'cluster_summary.tsv',sep='\t',index=False);pd.DataFrame(markers).to_csv(OUT/'within_stage_markers.tsv',sep='\t',index=False);pd.DataFrame(profiles).to_csv(OUT/'program_marker_profiles.tsv',sep='\t',index=False);pd.DataFrame(panelrows).to_csv(OUT/'program_panel_evidence.tsv',sep='\t',index=False);pd.DataFrame(mix).to_csv(OUT/'sample_mixing.tsv',sep='\t',index=False);pd.DataFrame(resolution).to_csv(OUT/'resolution_sensitivity.tsv',sep='\t',index=False);pd.DataFrame(globalstability).to_csv(OUT/'clustering_stability.tsv',sep='\t',index=False);pd.DataFrame(hvgrows).to_csv(OUT/'embedding_genes.tsv',sep='\t',index=False)
(OUT/'methods.json').write_text(json.dumps({'input':args.input,'population':'43249 previously ventricular-like nuclei','primary_resolution':.6,'resolution_checks':[.3,.6,.9,1.2],'seed':42,'seed_checks':[7,123],'HVG':'2000 sample-aware genes independently per age; mitochondrial/ribosomal and listed sex markers excluded','sex_embedding_exclusions':SEX,'PCA_components':30,'Harmony':'Within each age, sample key, theta2 lambda1 nclust30 max20iterations; no cross-age alignment','sensitivity':'Regress within-age QC/sex-marker covariates from Harmony PCs; can remove biological state and is not treated as proof of absence','markers':'Original corrected-count log1p/10000; cluster vs same-age rest, each library checked independently; no inferential tests','program_panels':PANELS,'support_heuristics':{'nuclei':100,'nuclei_per_library':30,'max_library_fraction':.9,'cross_library_markers':3,'seed_min_best_Jaccard':.5,'seed_mean_best_Jaccard':.65,'QC_sex_best_Jaccard':.5,'uncorrected_best_Jaccard':.4},'limitations':'Two library labels are not assumed to be independent animals. Higher resolution is not proof of new types. Anatomical program names remain candidates. Metadata-only H5ADs contain no expression. Count matrices unchanged.','versions':{k:importlib.metadata.version(k) for k in ['scanpy','anndata','harmonypy','numpy','scipy','scikit-learn']}},indent=2));print('DONE',flush=True)
