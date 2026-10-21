"""Determine whether Irx5/Kcnd2 associations depend on gate or existing state composition."""
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import spearmanr
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/ventricular_additional_axes'
f=ROOT/'results/ventricular_microenvironment/cell_annotations_and_scores.csv.gz';first=pd.read_csv(f,nrows=0).columns[0];cols=[first,'gate_origin','raw_total_counts','raw_n_genes_by_counts','decontx_contamination','raw_pct_counts_mt','cycling_marker_count'];meta=pd.read_csv(f,usecols=cols,index_col=0)
x=pd.read_csv(OUT/'targeted_log_expression.csv.gz',index_col=0);meta=meta.loc[x.index];rows=[]
def rho(a,b):return float(spearmanr(a,b).statistic) if np.std(a)>1e-10 and np.std(b)>1e-10 else 0.
for (age,sample),idx in x.groupby(['age','sample']).indices.items():
 for gate in ['expanded','original_core']:
  loc=idx if gate=='expanded' else idx[meta.iloc[idx].gate_origin.eq('core_ventricular').to_numpy()];d=x.iloc[loc];q=meta.iloc[loc];a=d.Irx5.to_numpy();b=d.Kcnd2.to_numpy();cov=np.column_stack([np.log1p(q.raw_total_counts),np.log1p(q.raw_n_genes_by_counts),q.decontx_contamination,q.raw_pct_counts_mt,q.cycling_marker_count]);cov=(cov-cov.mean(axis=0))/np.maximum(cov.std(axis=0),1e-8)
  def residual(v,c):return v-c@np.linalg.lstsq(c,v-v.mean(),rcond=None)[0]
  ac=residual(a,cov);bc=residual(b,cov);states=pd.get_dummies(d.microenvironment_cluster,dtype=float).to_numpy();states=states-states.mean(axis=0);combined=np.column_stack([cov,states]);asc=residual(a,combined);bsc=residual(b,combined)
  rows.append({'stage':age,'sample':sample,'gate':gate,'nuclei':len(d),'Irx5_Kcnd2_spearman':rho(a,b),'QC_cycle_residual_spearman':rho(ac,bc),'QC_cycle_and_existing_state_residual_spearman':rho(asc,bsc),'interpretation':'Conditioning on existing states can remove true spatial/electrical biological structure; not proof of absent within-state regulation.'})
pd.DataFrame(rows).to_csv(OUT/'electrical_axis_gate_state_sensitivity.tsv',sep='\t',index=False);print(pd.DataFrame(rows).round(3).to_string(index=False))
