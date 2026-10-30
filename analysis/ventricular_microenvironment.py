"""Unsupervised ventricular-enriched states interpreted against wall/oxygen/mechanics programs.
Programs are interpreted post hoc; primary unsupervised HVGs can include their genes. Orthogonal feature-exclusion, library and raw-count checks
are supplied. Expanded gate retains previously excluded Slit2-overlap reference CMs.
"""
from pathlib import Path
import os,argparse,json,re,importlib.metadata
os.environ.setdefault('MPLCONFIGDIR','/tmp/heart-mpl');os.environ.setdefault('NUMBA_CACHE_DIR','/tmp/heart-numba')
import numpy as np,pandas as pd,h5py,scanpy as sc,anndata as ad,harmonypy
from anndata.io import read_elem,sparse_dataset
from scipy import sparse
from scipy.stats import spearmanr
from sklearn.neighbors import NearestNeighbors
from sklearn.metrics import adjusted_rand_score
p=argparse.ArgumentParser();p.add_argument('--input',default='/workspace/heart-atlas/results/cell_identities_subtypes/snrna_qc_annotated_subclustered.h5ad');args=p.parse_args()
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/ventricular_microenvironment';OUT.mkdir(exist_ok=True);cfg=json.loads((OUT/'program_definitions.json').read_text());programs=cfg['programs'];sc.settings.n_jobs=4
SEX=['Xist','Tsix','Ddx3y','Kdm5d','Eif2s3y','Uty','Usp9y','Zfy1','Zfy2','Sry','Rbmy']
with h5py.File(args.input) as f:
 obs=read_elem(f['obs']);var=read_elem(f['var']);core=obs.proposed_subtype.eq('Ventricular-like');ref=obs.subcluster_id.isin(['CAR_3','CAR_9']);ix=np.flatnonzero((core|ref).to_numpy());obs=obs.iloc[ix].copy();obs['gate_origin']=np.where(core.iloc[ix],'core_ventricular','Slit2_overlap_reference')
 print('READ_CORRECTED',len(obs),flush=True);counts=sparse_dataset(f['X'])[ix,:].tocsr().astype(np.float32);print('READ_ORIGINAL',flush=True);raw=sparse_dataset(f['layers/counts_raw'])[ix,:].tocsr().astype(np.float32)
assert len(obs)==44463 and obs.gate_origin.eq('core_ventricular').sum()==43249 and obs.gate_origin.eq('Slit2_overlap_reference').sum()==1214
original_counts_totals=np.asarray(counts.sum(axis=1)).ravel();original_raw_totals=np.asarray(raw.sum(axis=1)).ravel()
def normalize(x):
 z=(sparse.diags(1e4/np.asarray(x.sum(axis=1)).ravel())@x).tocsr();z.data=np.log1p(z.data);return z
log=normalize(counts);rlog=normalize(raw);genes=var.gene_symbol.astype(str).to_numpy()
def mean(x,m):return np.asarray(x[m].mean(axis=0)).ravel()
def fraction(x,m):return np.asarray((x[m]>=1).mean(axis=0)).ravel()
def rho(x,y):return float(spearmanr(x,y).statistic) if np.std(x)>1e-10 and np.std(y)>1e-10 else 0.
def z_by_library(values,libs):
 out=np.empty_like(values,dtype=np.float64)
 for lib in np.unique(libs):
  m=libs==lib;out[m]=(values[m]-values[m].mean(axis=0))/np.maximum(values[m].std(axis=0),1e-8)
 return out
lookup={};mapping=[];all_program_features=[];feature_abundance=np.asarray(counts.sum(axis=0)).ravel()
for key,d in programs.items():
 js=[]
 for gene in d['genes']:
  found=np.flatnonzero(genes==gene);used=gene
  if not len(found):
   for alias in cfg.get('aliases',{}).get(gene,[]):
    found=np.flatnonzero(genes==alias)
    if len(found):used=alias;break
  all_program_features.extend(found)
  selected=found[np.argmax(feature_abundance[found]) : np.argmax(feature_abundance[found])+1] if len(found)>1 else found
  js.extend(selected);mapping.append({'program':key,'requested_gene':gene,'used_gene':used if len(found) else '', 'gene_ids':';'.join(var.index[found]),'score_gene_ids':';'.join(var.index[selected]),'feature_selection':'highest pooled corrected-count abundance if duplicated symbol','present':bool(len(found))})
 lookup[key]=np.unique(js)
allpanel=np.unique(all_program_features);panelmask=np.zeros(len(var),dtype=bool);panelmask[allpanel]=True
pd.DataFrame(mapping).to_csv(OUT/'program_gene_mapping.tsv',sep='\t',index=False)
prior=pd.read_csv(ROOT/'results/ventricular_within_stage/cell_annotations.csv.gz',index_col=0)
obs['prior_stage_group']=obs.index.map(prior.stage_subcluster_id.to_dict()).fillna('not_in_previous_ventricular_gate')
allann=[];summaries=[];gene_evidence=[];program_evidence=[];quality=[];smoothness=[];controls=[];metrics=[];resolutions=[];geneprofiles=[];stability=[];cohort=[];hvgrows=[];program_correlations=[]
def best_jaccard(ref,other,label):
 a=ref==label;matches=[(float(np.sum(a&(other==k))/np.sum(a|(other==k))),str(k)) for k in np.unique(other)];return max(matches,key=lambda x:x[0])
for stage in ['P0','P7','P14','P21']:
 print('STAGE',stage,flush=True);ind=np.flatnonzero(obs.age.eq(stage).to_numpy());o=obs.iloc[ind].copy()
 for col in o.select_dtypes('category').columns:o[col]=o[col].cat.remove_unused_categories()
 libs=o['sample'].astype(str).to_numpy();n=len(o);x=counts[ind].tocsr();rx=raw[ind].tocsr();l=log[ind].tocsr();rl=rlog[ind].tocsr();so=OUT/stage;so.mkdir(exist_ok=True)
 # Expression-matched controls are fixed from corrected-expression means and shared in raw sensitivity.
 avg=np.asarray(l.mean(axis=0)).ravel();valid=(avg>0)&~var.mt.to_numpy()&~var.ribo.to_numpy()&~var.gene_symbol.isin(SEX).to_numpy();valid_indices=np.flatnonzero(valid);ranks=np.argsort(np.argsort(avg[valid_indices],kind='stable'),kind='stable');bins=np.full(len(var),-1);bins[valid_indices]=np.minimum(19,ranks*20//len(valid_indices));score=[];rawscore=[];coverage=[];rawcoverage=[]
 for pi,(key,d) in enumerate(programs.items()):
  js=lookup[key];rng=np.random.default_rng(42+pi);chosen=[]
  for b in np.unique(bins[js]):
   if b<0:continue
   pool=np.flatnonzero((bins==b)&~panelmask&valid);chosen.extend(rng.choice(pool,min(50,len(pool)),replace=False))
  control=np.unique(chosen);assert len(js)>=2 and len(control)>0
  controls.append({'stage':stage,'program':key,'n_panel_features':len(js),'n_controls':len(control),'control_gene_ids':';'.join(var.index[control])})
  score.append(np.asarray(l[:,js].mean(axis=1)).ravel()-np.asarray(l[:,control].mean(axis=1)).ravel());rawscore.append(np.asarray(rl[:,js].mean(axis=1)).ravel()-np.asarray(rl[:,control].mean(axis=1)).ravel());coverage.append(np.asarray((x[:,js]>=1).mean(axis=1)).ravel());rawcoverage.append(np.asarray((rx[:,js]>=1).mean(axis=1)).ravel())
 score=np.array(score).T;rawscore=np.array(rawscore).T;coverage=np.array(coverage).T;rawcoverage=np.array(rawcoverage).T;z=z_by_library(score,libs);rz=z_by_library(rawscore,libs)
 covlist=[np.log1p(o.raw_total_counts.to_numpy()),np.log1p(o.raw_n_genes_by_counts.to_numpy()),o.decontx_contamination.to_numpy(),o.raw_pct_counts_mt.to_numpy(),o.cycling_marker_count.to_numpy()]
 for gene in ['Xist','Ddx3y','Kdm5d','Eif2s3y','Uty']:
  j=np.flatnonzero(genes==gene);covlist.append(np.asarray(l[:,j].sum(axis=1)).ravel())
 cov=np.array(covlist).T;cov=z_by_library(cov,libs);betas=np.linalg.lstsq(cov,z,rcond=None)[0];residual=z-cov@betas;adjusted_z=z_by_library(residual,libs)
 for pi,key in enumerate(programs):
  o['score_'+key]=score[:,pi];o['z_'+key]=z[:,pi];o['raw_z_'+key]=rz[:,pi];o['coverage_'+key]=coverage[:,pi];o['raw_coverage_'+key]=rawcoverage[:,pi];o['qc_adjusted_z_'+key]=adjusted_z[:,pi]
  quality.append({'stage':stage,'program':key,'features':len(lookup[key]),'median_detected_fraction':float(np.median(coverage[:,pi])),'fraction_cells_detecting_ge20pct_panel':float(np.mean(coverage[:,pi]>=.2)),'corrected_raw_score_spearman':rho(z[:,pi],rz[:,pi]),'QC_sex_cycle_explained_fraction':float(1-np.sum(residual[:,pi]**2)/np.sum(z[:,pi]**2)),'score_vs_log_raw_UMIs':rho(z[:,pi],cov[:,0]),'score_vs_log_raw_genes':rho(z[:,pi],cov[:,1]),'score_vs_contamination':rho(z[:,pi],cov[:,2]),'score_vs_mito':rho(z[:,pi],cov[:,3]),'score_vs_cycling':rho(z[:,pi],cov[:,4])})
 for i,key in enumerate(programs):
  for j,other in enumerate(programs):program_correlations.append({'stage':stage,'program_A':key,'program_B':other,'original_spearman':rho(z[:,i],z[:,j]),'QC_sex_cycle_residual_spearman':rho(adjusted_z[:,i],adjusted_z[:,j])})
 print('PCA_HARMONY',stage,flush=True)
 w=ad.AnnData(l.copy(),obs=o[['sample','age']].copy(),var=var.copy());sc.pp.highly_variable_genes(w,n_top_genes=2000,flavor='seurat',batch_key='sample');hv=w.var.highly_variable.to_numpy()&~var.mt.to_numpy()&~var.ribo.to_numpy()&~var.gene_symbol.isin(SEX).to_numpy()
 for j in np.flatnonzero(hv):hvgrows.append({'stage':stage,'gene_id':var.index[j],'gene_symbol':genes[j],'excluded_in_hypothesis_sensitivity':bool(panelmask[j])})
 representations={}
 for name,features in [('all_gene_HVG',hv),('hypothesis_genes_excluded',hv&~panelmask)]:
  t=w[:,features].copy();sc.pp.scale(t,zero_center=False,max_value=10);sc.pp.pca(t,n_comps=30,random_state=42,svd_solver='arpack');pc=t.obsm['X_pca'].copy();h=harmonypy.run_harmony(pc,o[['sample']],['sample'],theta=2,lamb=1,nclust=30,max_iter_harmony=20,max_iter_kmeans=10,random_state=42,ncores=4,verbose=False);cor=np.asarray(h.Z_corr);cor=cor if cor.shape==pc.shape else cor.T;cor=cor-cor.mean(axis=0)+pc.mean(axis=0);representations[name]=cor
  if name=='all_gene_HVG':representations['uncorrected']=pc.copy()
 representations['QC_sex_cycle_sensitivity']=representations['all_gene_HVG']-cov@np.linalg.lstsq(cov,representations['all_gene_HVG']-representations['all_gene_HVG'].mean(axis=0),rcond=None)[0]
 q=ad.AnnData(sparse.csr_matrix((n,0)),obs=o.copy());print('GRAPHS',stage,flush=True)
 for name,v in representations.items():
  q.obsm['X_'+name]=v;sc.pp.neighbors(q,n_neighbors=15,n_pcs=30,use_rep='X_'+name,random_state=42,key_added=name)
  for resolution in ([.3,.6,.9,1.2] if name=='all_gene_HVG' else [.6]):
   key=f'{name}_r{resolution:g}';sc.tl.leiden(q,neighbors_key=name,resolution=resolution,flavor='igraph',n_iterations=-1,random_state=42,key_added=key);resolutions.append({'stage':stage,'representation':name,'resolution':resolution,'clusters':int(q.obs[key].nunique())})
  nn=NearestNeighbors(n_neighbors=31,n_jobs=4).fit(v).kneighbors(v,return_distance=False)[:,1:];freq=o['sample'].value_counts().to_numpy()/n;expected=(1-(freq**2).sum())*n/(n-1);cross=float(np.mean(libs[:,None]!=libs[nn]));metrics.append({'stage':stage,'representation':name,'cross_library_neighbor_fraction':cross,'random_mixing_reference':float(expected),'observed_over_reference':float(cross/expected)})
 for seed in [7,123]:sc.tl.leiden(q,neighbors_key='all_gene_HVG',resolution=.6,flavor='igraph',n_iterations=-1,random_state=seed,key_added=f'seed_{seed}')
 sc.tl.umap(q,neighbors_key='all_gene_HVG',min_dist=.3,maxiter=250,random_state=42);q.obs['UMAP1']=q.obsm['X_umap'][:,0];q.obs['UMAP2']=q.obsm['X_umap'][:,1];q.obs['microenvironment_cluster']=stage+'_M'+q.obs['all_gene_HVG_r0.6'].astype(str)
 ref=q.obs['all_gene_HVG_r0.6'].astype(str).to_numpy();alternatives={key:q.obs[key].astype(str).to_numpy() for key in ['uncorrected_r0.6','hypothesis_genes_excluded_r0.6','QC_sex_cycle_sensitivity_r0.6','seed_7','seed_123','all_gene_HVG_r0.3','all_gene_HVG_r0.9','all_gene_HVG_r1.2']}
 for key,v in alternatives.items():stability.append({'stage':stage,'comparison':key,'adjusted_rand_index':adjusted_rand_score(ref,v)})
 core=q.obs.gate_origin.eq('core_ventricular');stability.append({'stage':stage,'comparison':'core_gate_vs_previous_core_only_r0.6','adjusted_rand_index':adjusted_rand_score(q.obs.loc[core,'microenvironment_cluster'],q.obs.loc[core,'prior_stage_group'])})
 # Neighborhood smoothness versus within-library shuffled references, not animal-level tests.
 graph=q.obsp['all_gene_HVG_connectivities'].tocsr();rowsum=np.asarray(graph.sum(axis=1)).ravel();neighbor=(graph@z)/np.maximum(rowsum[:,None],1e-8);neighbor_adj=(graph@adjusted_z)/np.maximum(rowsum[:,None],1e-8)
 nulls=[[] for _ in programs];rng=np.random.default_rng(42)
 for trial in range(20):
  perm=np.arange(n)
  for lib in np.unique(libs):
   loc=np.flatnonzero(libs==lib);perm[loc]=rng.permutation(loc)
  shuffled=z[perm];near=(graph@shuffled)/np.maximum(rowsum[:,None],1e-8)
  for pi in range(len(programs)):nulls[pi].append(rho(shuffled[:,pi],near[:,pi]))
 excluded_graph=q.obsp['hypothesis_genes_excluded_connectivities'].tocsr();excluded_sum=np.asarray(excluded_graph.sum(axis=1)).ravel();excluded_neighbor=(excluded_graph@z)/np.maximum(excluded_sum[:,None],1e-8)
 for pi,key in enumerate(programs):smoothness.append({'stage':stage,'program':key,'score_neighbor_spearman':rho(z[:,pi],neighbor[:,pi]),'QC_sex_cycle_residual_neighbor_spearman':rho(adjusted_z[:,pi],neighbor_adj[:,pi]),'interpreted_genes_excluded_neighbor_spearman':rho(z[:,pi],excluded_neighbor[:,pi]),'within_library_shuffle_reference_mean':float(np.mean(nulls[pi])),'shuffle_reference_max':float(np.max(nulls[pi])),'interpretation':'Neighborhood coherence, not a spatial measurement or proof of discrete types'})
 # All expression evidence remains on original corrected/original integer counts, never PCs.
 clusterids=sorted(q.obs.microenvironment_cluster.unique(),key=lambda s:int(s.split('M')[-1]));codes=pd.Categorical(q.obs.microenvironment_cluster,categories=clusterids).codes;g=sparse.csr_matrix((np.ones(n,dtype=np.float64),(codes,np.arange(n))),shape=(len(clusterids),n));ns=np.asarray(g.sum(axis=1)).ravel();logs=(g@l).toarray();rawlogs=(g@rl).toarray();fr=(g@(x>=1)).toarray()/ns[:,None];rfr=(g@(rx>=1)).toarray()/ns[:,None];sums=(g@x).toarray();libdata={}
 for lib in np.unique(libs):
  lm=libs==lib;libdata[lib]=(np.asarray(g[:,lm].sum(axis=1)).ravel(),(g[:,lm]@l[lm]).toarray(),(g[:,lm]@rl[lm]).toarray())
 for ci,cid in enumerate(clusterids):
  m=q.obs.microenvironment_cluster.eq(cid).to_numpy();b=~m;sub=q.obs.loc[m];dm=logs[ci]/ns[ci]-(logs.sum(axis=0)-logs[ci])/(n-ns[ci]);rdm=rawlogs[ci]/ns[ci]-(rawlogs.sum(axis=0)-rawlogs[ci])/(n-ns[ci]);contrasts={};rawcontrasts={}
  for lib,(lns,ls,rs) in libdata.items():
   a=lns[ci];bb=lns.sum()-a
   if min(a,bb)>=30:contrasts[lib]=ls[ci]/a-(ls.sum(axis=0)-ls[ci])/bb;rawcontrasts[lib]=rs[ci]/a-(rs.sum(axis=0)-rs[ci])/bb
  gene_support=np.sum(np.array([((contrasts[k]>.05)&(rawcontrasts[k]>.05)) for k in contrasts]),axis=0) if contrasts else np.zeros(len(var),dtype=int)
  restfr=(fr*ns[:,None]).sum(axis=0)-fr[ci]*ns[ci];restfr/=n-ns[ci]
  eligible=(dm>.15)&(fr[ci]>=.15)&((fr[ci]-restfr)>=.05)&~var.mt.to_numpy()&~var.ribo.to_numpy()&~var.gene_symbol.isin(SEX).to_numpy();order=np.flatnonzero(eligible);order=order[np.argsort(dm[order])[::-1]][:80];supported=[];seen_symbols=set()
  for j in order:
   if gene_support[j]>=2 and genes[j] not in seen_symbols:supported.append(j);seen_symbols.add(genes[j])
  for rank,j in enumerate(order,1):gene_evidence.append({'stage':stage,'cluster':cid,'rank':rank,'gene_id':var.index[j],'gene_symbol':genes[j],'corrected_log_difference':dm[j],'raw_log_difference':rdm[j],'corrected_detected_fraction':fr[ci,j],'raw_detected_fraction':rfr[ci,j],'eligible_libraries':len(contrasts),'libraries_with_positive_corrected_and_raw_effect':int(gene_support[j])})
  candidate_programs=[];program_lookup={}
  for pi,(key,d) in enumerate(programs.items()):
   js=lookup[key];positive=[j for j in js if dm[j]>.1 and fr[ci,j]>=.1 and gene_support[j]>=2];anchors=[gene for gene in d['anchors'] if any((genes[j]==gene) for j in positive)];scoreeffect=float(z[m,pi].mean()-z[b,pi].mean());raweffect=float(rz[m,pi].mean()-rz[b,pi].mean());qceffect=float(adjusted_z[m,pi].mean()-adjusted_z[b,pi].mean());perlib=[];rawperlib=[];qcperlib=[]
   for lib in contrasts:
    ma=m&(libs==lib);mb=b&(libs==lib);perlib.append(float(z[ma,pi].mean()-z[mb,pi].mean()));rawperlib.append(float(rz[ma,pi].mean()-rz[mb,pi].mean()));qcperlib.append(float(adjusted_z[ma,pi].mean()-adjusted_z[mb,pi].mean()))
   passgene=len(positive)>=d['min_supported_genes'] and len(anchors)>=d['anchor_min'];passed=passgene and scoreeffect>.35 and len(perlib)==2 and min(perlib)>.15 and min(rawperlib)>0 and raweffect>0
   row={'stage':stage,'cluster':cid,'program':key,'program_label':d['display_label'],'family':d['family'],'corrected_z_effect':scoreeffect,'raw_z_effect':raweffect,'QC_adjusted_z_effect':qceffect,'mean_corrected_detection_fraction':float(coverage[m,pi].mean()),'mean_raw_detection_fraction':float(rawcoverage[m,pi].mean()),'supported_gene_count':len(positive),'supported_genes':';'.join(genes[positive]),'supported_anchors':';'.join(anchors),'library_corrected_effects':';'.join(map(str,perlib)),'library_raw_effects':';'.join(map(str,rawperlib)),'library_QC_adjusted_effects':';'.join(map(str,qcperlib)),'program_rule_pass':bool(passed),'QC_adjusted_program_support':bool(len(qcperlib)==2 and min(qcperlib)>.1)};program_evidence.append(row);program_lookup[key]=row
   if passed and d['family'] in ['hypothesis','metabolic']:candidate_programs.append((key,scoreeffect))
   for j in js:geneprofiles.append({'stage':stage,'cluster':cid,'program':key,'gene':genes[j],'mean_log1p_10k':logs[ci,j]/ns[ci],'corrected_detected_fraction':fr[ci,j],'raw_detected_fraction':rfr[ci,j],'corrected_log_difference':dm[j],'raw_log_difference':rdm[j],'positive_corrected_raw_libraries':int(gene_support[j])})
  # A HIF label additionally requires the non-glycolytic sensitivity, not glycolysis alone.
  candidate_programs=[v for v in candidate_programs if v[0]!='HIF_response' or program_lookup['HIF_non_glycolytic']['program_rule_pass']]
  candidate_programs.sort(key=lambda v:v[1],reverse=True)
  if program_lookup['Cycling']['program_rule_pass'] and sub.cycling_state.eq('cycling-marker positive').mean()>=.1:best='Cycling';label='Cycling-associated control state'
  elif candidate_programs:best=candidate_programs[0][0];label=programs[best]['display_label']
  else:
   best='';names=[genes[j] for j in supported if not re.match(r'^(Gm\d|\d)',genes[j])];label='Marker-defined: '+ '/'.join(names[:2]) if names else 'Uncharacterized ventricular-enriched group'
  rlabel=str(q.obs.loc[m,'all_gene_HVG_r0.6'].iloc[0]);matches={k:best_jaccard(ref,v,rlabel) for k,v in alternatives.items()};seeds=[matches['seed_7'][0],matches['seed_123'][0]];library_counts=sub['sample'].value_counts();sizepass=len(sub)>=100 and (library_counts>=30).sum()==2 and library_counts.max()/len(sub)<=.9;seedpass=min(seeds)>=.5 and np.mean(seeds)>=.65;markerpass=len(supported)>=3;replicate=sizepass and seedpass and markerpass
  nuisance=matches['QC_sex_cycle_sensitivity_r0.6'][0]<.5;paneldependent=matches['hypothesis_genes_excluded_r0.6'][0]<.5;sampledependent=matches['uncorrected_r0.6'][0]<.4
  status='supported RNA state' if replicate and not(nuisance or paneldependent or sampledependent) else 'replicate-supported; sensitivity-dependent' if replicate else 'candidate/review'
  flags=[]
  for condition,name in [(not sizepass,'size/library'),(not seedpass,'seed-sensitive'),(not markerpass,'limited cross-library raw/corrected markers'),(nuisance,'QC/sex/cycling-sensitive'),(paneldependent,'interpreted-gene-dependent'),(sampledependent,'sample-adjustment-sensitive')]:
   if condition:flags.append(name)
  summaries.append({'stage':stage,'microenvironment_cluster':cid,'nuclei':len(sub),'microenvironment_label':label,'selected_program':best,'support_status':status,'review_flags':';'.join(flags),'core_ventricular_n':int(sub.gate_origin.eq('core_ventricular').sum()),'Slit2_overlap_reference_n':int(sub.gate_origin.eq('Slit2_overlap_reference').sum()),'libraries_ge30':int((library_counts>=30).sum()),'largest_library_fraction':float(library_counts.max()/len(sub)),'both_library_raw_corrected_marker_count':len(supported),'both_library_markers':';'.join(genes[supported[:15]]),'seed_mean_best_jaccard':float(np.mean(seeds)),'seed_min_best_jaccard':float(min(seeds)),'QC_sensitivity_best_jaccard':matches['QC_sex_cycle_sensitivity_r0.6'][0],'panel_excluded_best_jaccard':matches['hypothesis_genes_excluded_r0.6'][0],'uncorrected_best_jaccard':matches['uncorrected_r0.6'][0],'cycling_marker_positive_fraction':float(sub.cycling_state.eq('cycling-marker positive').mean()),'median_raw_UMIs':float(sub.raw_total_counts.median()),'median_contamination':float(sub.decontx_contamination.median())})
  q.obs.loc[m,'microenvironment_label']=label;q.obs.loc[m,'support_status']=status
 for col in ['gate_origin','sample','subcluster_id','prior_stage_group']:pd.crosstab(q.obs.microenvironment_cluster,q.obs[col]).to_csv(so/f'cluster_by_{col}.tsv',sep='\t')
 pb=pd.DataFrame(sums.T,index=var.index,columns=clusterids);pb.insert(0,'gene_symbol',genes);pb.index.name='gene_id';pb.to_csv(so/'corrected_pseudobulk.tsv.gz',sep='\t',float_format='%.6g')
 q.obs.to_csv(so/'cell_annotations_and_scores.csv.gz');q.uns['analysis_method']={'gate':'Expanded ventricular-enriched: prior ventricular-like plus ambiguous Slit2/conduction-overlap references; prior atrial-like excluded.','discovery':'Unsupervised within-stage HVG/PCA/Harmony; programs used after clustering.','coordinates':'Independent per stage; not comparable across ages.','gene_expression':'Absent in metadata-only H5AD. Original corrected/raw counts unaltered.','scores':'Expression-matched control genes; within-library standardization; raw and QC/sex/cycle sensitivities. No direct spatial/oxygen/mechanical measurements.'};q.write_h5ad(so/'embeddings_and_scores_metadata_only.h5ad',compression='gzip');allann.append(q.obs);cohort.append({'stage':stage,'expanded_nuclei':n,'core_nuclei':int(o.gate_origin.eq('core_ventricular').sum()),'reference_nuclei':int(o.gate_origin.eq('Slit2_overlap_reference').sum()),'primary_clusters':len(clusterids),'corrected_count_sum':float(x.sum()),'original_integer_count_sum':float(rx.sum())});print('FINISHED',stage,len(clusterids),flush=True)
for name,data in [('cluster_summary',summaries),('cluster_gene_markers',gene_evidence),('cluster_program_evidence',program_evidence),('module_score_quality',quality),('gradient_metrics',smoothness),('matched_control_genes',controls),('sample_mixing',metrics),('resolution_sensitivity',resolutions),('program_gene_profiles',geneprofiles),('clustering_stability',stability),('cohort_summary',cohort),('embedding_genes',hvgrows),('program_correlations',program_correlations)]:pd.DataFrame(data).to_csv(OUT/(name+'.tsv'),sep='\t',index=False)
pd.concat(allann).to_csv(OUT/'cell_annotations_and_scores.csv.gz')
(OUT/'methods.json').write_text(json.dumps({'input':args.input,'expanded_nuclei':44463,'core_nuclei':43249,'reference_nuclei':1214,'primary_resolution':.6,'resolution_checks':[.3,.6,.9,1.2],'PCA':30,'neighbors':15,'UMAP_min_dist':.3,'seed':42,'seed_checks':[7,123],'scoring':'Per-stage mean log-normalized program minus matched-expression control means (20 gene bins; up to 50 controls/bin; exclude all interpreted program genes from controls); control genes shared for raw-count scoring; scores standardized separately in each original library.','raw_sensitivity':'Original integer counts normalized independently with original library totals; no rounding of corrected fractions.','duplicate_gene_symbols':'Program scores select one feature per symbol, choosing highest pooled corrected-count abundance; support counts distinct symbols. All same-symbol features excluded in panel sensitivity; original counts remain unchanged.','gene_exclusion_sensitivity':'Remove every interpreted program feature from the stage HVG embedding; repeat PCA/Harmony/graph/clustering. Not a proof of absent biology if neighborhoods change.','QC_score_sensitivity':'Within-library-centered log raw UMIs, log raw detected genes, contamination, mitochondrial percentage, cycling-marker count and Xist/Y-marker expression; can remove biological state.','program_labeling':'Prespecified source/marker conditions, >=.35 z effect, >.15 in both original libraries, raw effects positive, >=3 or panel-specific positive genes with corrected and raw positive contrasts in both libraries; HIF also requires non-glycolytic program.','state_support':'n>=100, >=30 each library, largest fraction<=.9, >=3 both-library corrected/raw markers; mean seed Jaccard>=.65, min>=.5; QC and feature-exclusion Jaccard>=.5, uncorrected Jaccard>=.4. These are descriptive heuristics.','gradient_check':'Score-neighbor Spearman on expression graph and after QC adjustment; 20 within-library shuffled reference distributions, not animal-level inference or spatial distance.','gate_correction':'Earlier Etv1/Gja5-overlap label could exclude paper Slit2 trabecular CMs; retained as ambiguous reference, not assigned Purkinje identity.','limits':'No measured wall position, tissue oxygen, force, metabolic flux, independent animal replication or causal environmental attribution. Author Table3 panels from same dataset are not independent validation.','versions':{k:importlib.metadata.version(k) for k in ['scanpy','anndata','harmonypy','numpy','scipy','scikit-learn']}},indent=2));print('DONE',flush=True)
