"""Reconcile conservative program labels without changing embeddings or clusters."""
from pathlib import Path
import re
import pandas as pd,anndata as ad
ROOT=Path(__file__).resolve().parents[1];out=ROOT/'results/ventricular_within_stage'
s=pd.read_csv(out/'cluster_summary.tsv',sep='\t');p=pd.read_csv(out/'program_panel_evidence.tsv',sep='\t');m=pd.read_csv(out/'within_stage_markers.tsv',sep='\t')
p['panel_rule_pass']=p.panel_rule_pass&(p.relative_score>0)
for i,row in s.iterrows():
 panel=p[p.stage_subcluster_id.eq(row.stage_subcluster_id)];passed=panel[panel.panel_rule_pass];cycle=passed[passed.program.eq('Cycling/DNA-replication program')]
 others=passed[~passed.program.isin(['Cycling/DNA-replication program','Conduction-overlap check'])].sort_values('relative_score',ascending=False)
 if row.cycling_positive_fraction>=.1 and len(cycle):label='Cycling/DNA-replication program'
 elif len(others):label=others.program.iloc[0]
 else:
  g=m[m.stage_subcluster_id.eq(row.stage_subcluster_id)].sort_values('rank');valid=lambda xs:[gene for gene in xs if not re.match(r'^(Gm\d|\d)|Rik$',gene)]
  names=valid(g.loc[g.positive_libraries.ge(2),'gene_symbol']) or valid(g.gene_symbol)
  label='Marker-defined: '+('/'.join(names[:2])) if names else 'Uncharacterized ventricular state'
 selected_program=label if label in set(panel.program) else ''
 selected_genes=set(str(panel.loc[panel.program.eq(label),'supported_panel_genes'].iloc[0]).split(';')) if selected_program else set()
 replicated=set(m.loc[m.stage_subcluster_id.eq(row.stage_subcluster_id)&m.positive_libraries.ge(2),'gene_symbol'])
 if label=='Calcium/contractile program' and selected_genes=={'Srl','Actn2'}:label='Srl/Actn2-enriched structural program'
 elif label=='Adhesion/structural program' and len(selected_genes)<=2 and {'Fgf14','Thsd4'}.issubset(replicated):label='Fgf14/Thsd4-enriched state'
 elif label=='Neonatal contractile program' and {'Myh7','Nppa'}.issubset(selected_genes):label='Myh7/Nppa-enriched program'
 s.loc[i,'selected_program_panel']=selected_program
 status='supported transcriptional state' if row.replicate_program_support and not row.qc_sex_sensitive and not row.sample_adjustment_sensitive else 'replicate-supported; adjustment-sensitive' if row.replicate_program_support else 'candidate/review'
 if 'candidate' in label.lower():status='candidate anatomical association; '+status
 s.loc[i,'stage_state_label']=label;s.loc[i,'stage_support_status']=status
mapping=s.set_index('stage_subcluster_id');ann=pd.read_csv(out/'cell_annotations.csv.gz',index_col=0)
for col in ['stage_state_label','stage_support_status']:ann[col]=ann.stage_subcluster_id.map(mapping[col])
s.to_csv(out/'cluster_summary.tsv',sep='\t',index=False);p.to_csv(out/'program_panel_evidence.tsv',sep='\t',index=False);ann.to_csv(out/'cell_annotations.csv.gz')
for stage in ['P0','P7','P14','P21']:
 q=ad.read_h5ad(out/stage/'embeddings_metadata_only.h5ad')
 for col in ['stage_state_label','stage_support_status']:q.obs[col]=q.obs.stage_subcluster_id.astype(str).map(mapping[col]).astype(str)
 q.write_h5ad(out/stage/'embeddings_metadata_only.h5ad',compression='gzip');q.obs.to_csv(out/stage/'cell_annotations.csv.gz')
print(s[['stage_subcluster_id','nuclei','stage_state_label','stage_support_status']].to_string(index=False))
