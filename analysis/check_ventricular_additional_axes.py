"""Literature-guided gene checks without imposing functional/anatomical cell labels."""
from pathlib import Path
import json,argparse
import numpy as np,pandas as pd,h5py
from anndata.io import read_elem,sparse_dataset
from scipy.stats import spearmanr
p=argparse.ArgumentParser();p.add_argument('--input',default='/workspace/heart-atlas/results/cell_identities_subtypes/snrna_qc_annotated_subclustered.h5ad');args=p.parse_args()
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/ventricular_additional_axes';OUT.mkdir(exist_ok=True)
GENES='Irx5 Irx3 Irx4 Kcnd2 Kcnd3 Kcnip2 Smyd1 Gja1 Gja5 Scn5a Scn10a Kcnj2 Kcnj12 Kcnb1 Kcna4 Kcna5 Kcnq1 Kcne1 Kcnh2 Hand1 Hand2 Tbx5 Isl1 Mfap2 Mcam Prelid2 Cdh13 Dlc1 Ebf2 Sox5 Ebf1 Kcnab1 Fhl1 Dusp27 Xirp1 Xirp2 Cnn1 Nexn Ndufb11 Ndufa4 Cox7c Cox5b Cryab Myh6 Myh7 Myl2 Tnnt2 Tnni3 Tnni1 Ryr2 Atp2a2 Casq2 Jph2 Srl Pln Cacna1c Cacna1d Tnni3k Ect2 Anln Aurkb Mki67 Top2a Fabp3 Cox6c Oxct1 Myl7 Hmgn2 Pgam1 Arntl Clock Per1 Per2 Per3 Cry1 Cry2 Nr1d1 Nfil3 Slit2 Ankrd1 Ddc Etv1 Cntn2 Nppa Nppb Hey2 Tbx20 Mycn Col1a1 Col1a2 Col3a1 Dcn Lum Fn1'.split()
ann=pd.read_csv(ROOT/'results/ventricular_microenvironment/cell_annotations_and_scores.csv.gz',index_col=0)
with h5py.File(args.input) as f:
 allobs=read_elem(f['obs']);var=read_elem(f['var']);rows=allobs.index.get_indexer(ann.index);assert (rows>=0).all();names=var.gene_symbol.astype(str).to_numpy()
 print('READ_TARGETED_CORRECTED',len(ann),flush=True);x=sparse_dataset(f['X'])[rows,:].tocsr().astype(np.float32);abundance=np.asarray(x.sum(axis=0)).ravel();selection=[];mapping=[]
 for gene in GENES:
  found=np.flatnonzero(names==gene)
  if len(found):
   j=int(found[np.argmax(abundance[found])]);selection.append(j);mapping.append({'gene':gene,'selected_gene_id':var.index[j],'all_matching_gene_ids':';'.join(var.index[found]),'present':True})
  else:mapping.append({'gene':gene,'selected_gene_id':'','all_matching_gene_ids':'','present':False})
 symbols=names[selection];corrected=x[:,selection].toarray();totals=np.asarray(x.sum(axis=1)).ravel();del x
 print('READ_TARGETED_RAW',flush=True);r=sparse_dataset(f['layers/counts_raw'])[rows,:].tocsr();raw=r[:,selection].toarray();rawtotals=np.asarray(r.sum(axis=1)).ravel();del r
log=np.log1p(corrected*1e4/totals[:,None]);rawlog=np.log1p(raw*1e4/rawtotals[:,None]);det=corrected>=1;rdet=raw>=1
profiles=[]
for fields in [['age'],['age','sample'],['age','microenvironment_cluster'],['age','sample','microenvironment_cluster']]:
 groups=ann.groupby(fields,observed=True).indices
 for key,idx in groups.items():
  key=key if isinstance(key,tuple) else (key,);meta=dict(zip(fields,key));n=len(idx)
  for j,gene in enumerate(symbols):profiles.append({**meta,'grouping':'+'.join(fields),'nuclei':n,'gene':gene,'mean_log1p_10k':float(log[idx,j].mean()),'mean_raw_log1p_10k':float(rawlog[idx,j].mean()),'fraction_corrected_ge1':float(det[idx,j].mean()),'fraction_raw_ge1':float(rdet[idx,j].mean())})
pd.DataFrame(profiles).to_csv(OUT/'gene_profiles_by_stage_library_cluster.tsv',sep='\t',index=False)
pd.DataFrame(mapping).to_csv(OUT/'gene_mapping.tsv',sep='\t',index=False)
lookup={gene:j for j,gene in enumerate(symbols)}
def rho(a,b):return float(spearmanr(a,b).statistic) if np.std(a)>1e-10 and np.std(b)>1e-10 else 0.
correlations=[]
for (stage,lib),idx in ann.groupby(['age','sample']).indices.items():
 sub=ann.iloc[idx];cov=np.column_stack([np.log1p(sub.raw_total_counts),np.log1p(sub.raw_n_genes_by_counts),sub.decontx_contamination,sub.raw_pct_counts_mt,sub.cycling_marker_count]);cov=(cov-cov.mean(axis=0))/np.maximum(cov.std(axis=0),1e-8)
 for g1,g2 in [('Irx5','Kcnd2'),('Kcnd2','Kcnip2'),('Kcnd2','Tnni3'),('Irx5','Smyd1'),('Hand1','Tbx5'),('Mcam','Mfap2')]:
  a=log[idx,lookup[g1]];b=log[idx,lookup[g2]];ar=a-cov@np.linalg.lstsq(cov,a-a.mean(),rcond=None)[0];br=b-cov@np.linalg.lstsq(cov,b-b.mean(),rcond=None)[0]
  # Cluster-mean association is descriptive and based on few groups; not biological replication.
  groups=sub.groupby('microenvironment_cluster').indices;ma=np.array([a[j].mean() for j in groups.values() if len(j)>=30]);mb=np.array([b[j].mean() for j in groups.values() if len(j)>=30])
  correlations.append({'stage':stage,'sample':lib,'gene_A':g1,'gene_B':g2,'nuclei':len(idx),'nucleus_spearman':rho(a,b),'QC_cycle_residual_spearman':rho(ar,br),'raw_nucleus_spearman':rho(rawlog[idx,lookup[g1]],rawlog[idx,lookup[g2]]),'clusters_with_ge30_nuclei':len(ma),'cluster_mean_spearman':rho(ma,mb),'joint_corrected_detection_fraction':float((det[idx,lookup[g1]]&det[idx,lookup[g2]]).mean())})
pd.DataFrame(correlations).to_csv(OUT/'within_library_gene_associations.tsv',sep='\t',index=False)
cm=['Tnnt2','Tnni3','Myh6','Myl2'];ecm=['Col1a1','Col3a1','Dcn'];cm_n=det[:,[lookup[g] for g in cm]].sum(axis=1);ecm_n=det[:,[lookup[g] for g in ecm]].sum(axis=1);raw_ecm=rdet[:,[lookup[g] for g in ecm]].sum(axis=1)
flags=pd.DataFrame({'age':ann.age,'sample':ann['sample'],'microenvironment_cluster':ann.microenvironment_cluster,'CM_markers_detected':cm_n,'ECM_markers_detected_corrected':ecm_n,'ECM_markers_detected_raw':raw_ecm,'ECM2_and_CM2_corrected':(ecm_n>=2)&(cm_n>=2),'ECM3_and_CM2_corrected':(ecm_n>=3)&(cm_n>=2),'decontx_contamination':ann.decontx_contamination,'scdblfinder_score':ann.scdblfinder_score},index=ann.index)
flags.groupby(['age','sample']).agg(nuclei=('age','size'),ECM2_CM2=('ECM2_and_CM2_corrected','sum'),ECM3_CM2=('ECM3_and_CM2_corrected','sum'),raw_ECM2=('ECM_markers_detected_raw',lambda v:int((v>=2).sum()))).to_csv(OUT/'ECM_CM_coexpression_counts.tsv',sep='\t')
flags[flags.ECM2_and_CM2_corrected].to_csv(OUT/'ECM_coexpression_candidates.csv.gz')
# Small per-cell expression table supports independent inspection; counts are original values, unadjusted.
expression=pd.DataFrame(log,index=ann.index,columns=symbols);expression.insert(0,'sample',ann['sample']);expression.insert(0,'age',ann.age);expression.insert(0,'microenvironment_cluster',ann.microenvironment_cluster);expression.to_csv(OUT/'targeted_log_expression.csv.gz',float_format='%.6g')
(OUT/'methods.json').write_text(json.dumps({'input':args.input,'nuclei':len(ann),'genes_requested':GENES,'gene_selection':'One feature per symbol; greatest pooled corrected abundance if duplicates. All original matrices unchanged.','expression':'log1p(count/whole-nucleus count total*10000), separately for corrected fractional and original integer counts. Detection>=1 count.','comparison':'Within age and original library; descriptive associations, no nucleus-level p-values or independent-animal inference.','residual':'Within-library linear removal of log RNA depth, log detected genes, contamination, mitochondrial fraction and cycling-marker count; conditional sensitivity, may remove biology.','ECM':'Co-detection is candidate evidence only; ambient RNA and remaining doublets remain alternatives. Does not establish ECM-producing CM subtype.','limits':'RNA is not ion current/action potential, wall position, chamber identity, DNA content/nucleation, clock phase or causal mechanism. Cross-species/developmental marker transfer requires validation.'},indent=2));print('DONE',flush=True)
