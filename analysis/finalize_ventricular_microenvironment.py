"""Add core-gate/orthogonal-gradient diagnostics and explicit cycling-control labels."""
from pathlib import Path
import json
import numpy as np,pandas as pd,anndata as ad
from scipy.stats import spearmanr
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/ventricular_microenvironment';cfg=json.loads((OUT/'program_definitions.json').read_text());programs=cfg['programs']
ann=pd.read_csv(OUT/'cell_annotations_and_scores.csv.gz',index_col=0);summary=pd.read_csv(OUT/'cluster_summary.tsv',sep='\t');evidence=pd.read_csv(OUT/'cluster_program_evidence.tsv',sep='\t');gradient=pd.read_csv(OUT/'gradient_metrics.tsv',sep='\t');gate=[]
def rho(a,b):return float(spearmanr(a,b).statistic) if np.std(a)>1e-10 and np.std(b)>1e-10 else 0.
def standardized(v,libs):
 out=np.zeros_like(v,dtype=float)
 for lib in np.unique(libs):
  m=libs==lib;out[m]=(v[m]-v[m].mean())/max(v[m].std(),1e-8)
 return out
for i,row in summary.iterrows():
 cycle=evidence[(evidence.cluster==row.microenvironment_cluster)&evidence.program.eq('Cycling')]
 if bool(cycle.program_rule_pass.iloc[0]) and row.cycling_marker_positive_fraction>=.1:
  summary.loc[i,'selected_program']='Cycling';summary.loc[i,'microenvironment_label']='Cycling-associated control state'
labelmap=summary.set_index('microenvironment_cluster').microenvironment_label
ann['microenvironment_label']=ann.microenvironment_cluster.map(labelmap)
for stage in ['P0','P7','P14','P21']:
 q=ad.read_h5ad(OUT/stage/'embeddings_and_scores_metadata_only.h5ad');q.obs['microenvironment_label']=q.obs.microenvironment_cluster.astype(str).map(labelmap).astype(str)
 old=ad.read_h5ad(ROOT/'results/ventricular_within_stage'/stage/'embeddings_metadata_only.h5ad');core=q.obs.loc[old.obs_names];libs=core['sample'].astype(str).to_numpy();g=old.obsp['harmony_connectivities'].tocsr();den=np.asarray(g.sum(axis=1)).ravel();panel_graph=q.obsp['hypothesis_genes_excluded_connectivities'].tocsr();pden=np.asarray(panel_graph.sum(axis=1)).ravel();ref=q.obs.gate_origin.eq('Slit2_overlap_reference')
 for key in programs:
  cs=standardized(core['score_'+key].to_numpy(),libs);cq=standardized(core['qc_adjusted_z_'+key].to_numpy(),libs)
  expanded=q.obs['z_'+key].to_numpy();resid=q.obs['qc_adjusted_z_'+key].to_numpy();m=(gradient.stage==stage)&gradient.program.eq(key)
  gradient.loc[m,'interpreted_genes_excluded_neighbor_spearman']=rho(expanded,(panel_graph@expanded)/np.maximum(pden,1e-8));gradient.loc[m,'interpreted_genes_excluded_residual_neighbor_spearman']=rho(resid,(panel_graph@resid)/np.maximum(pden,1e-8))
  gate.append({'stage':stage,'program':key,'core_nuclei':len(core),'reference_nuclei':int(ref.sum()),'reference_minus_core_expanded_z':float(q.obs.loc[ref,'z_'+key].mean()-q.obs.loc[~ref,'z_'+key].mean()),'core_only_graph_neighbor_spearman':rho(cs,(g@cs)/np.maximum(den,1e-8)),'core_only_graph_QC_residual_neighbor_spearman':rho(cq,(g@cq)/np.maximum(den,1e-8)),'interpretation':'Original core-only gene/Harmony graph; common expanded-gate controls, standardized again within each core library. Core residual uses expanded-gate QC fit.'})
 q.write_h5ad(OUT/stage/'embeddings_and_scores_metadata_only.h5ad',compression='gzip');q.obs.to_csv(OUT/stage/'cell_annotations_and_scores.csv.gz')
summary.to_csv(OUT/'cluster_summary.tsv',sep='\t',index=False);ann.to_csv(OUT/'cell_annotations_and_scores.csv.gz');gradient.to_csv(OUT/'gradient_metrics.tsv',sep='\t',index=False);pd.DataFrame(gate).to_csv(OUT/'core_gate_sensitivity.tsv',sep='\t',index=False)
methods=json.loads((OUT/'methods.json').read_text());methods['core_gate_sensitivity']='Reuse original 43,249-cell stage-specific core-only Harmony graphs; common program controls and core-library standardization; compare neighborhood coherence with expanded-gate results. Expanded-gate QC fit used for residual comparison.';methods['primary_feature_clarification']='Primary unsupervised HVGs can contain interpreted genes; programs are not used as a handpicked-only feature basis. Separate all-program-gene-excluded embedding tests dependence.';methods['cycling_labels']='Multi-gene cycling-enriched clusters explicitly labeled as control states, not microenvironment-induced types.';(OUT/'methods.json').write_text(json.dumps(methods,indent=2));print('Finalized conservative labels, core-gate comparisons and gene-excluded neighborhood checks.')
