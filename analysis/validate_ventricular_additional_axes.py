from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/ventricular_additional_axes'
e=pd.read_csv(OUT/'targeted_log_expression.csv.gz',index_col=0);source=pd.read_csv(ROOT/'results/ventricular_microenvironment/cell_annotations_and_scores.csv.gz',usecols=lambda c:c in ['cell_id','age','sample','microenvironment_cluster'])
# Original index name is retained; use the source annotation index for an unambiguous population check.
source_index=pd.read_csv(ROOT/'results/ventricular_microenvironment/cell_annotations_and_scores.csv.gz',usecols=[0],index_col=0).index
assert len(e)==44463 and e.index.is_unique and e.index.equals(source_index)
map=pd.read_csv(OUT/'gene_mapping.tsv',sep='\t');assert map.present.all() and map.gene.is_unique;assert set(e.columns[3:])==set(map.gene)
assert np.isfinite(e.iloc[:,3:].to_numpy()).all() and (e.iloc[:,3:].to_numpy()>=0).all()
p=pd.read_csv(OUT/'gene_profiles_by_stage_library_cluster.tsv',sep='\t');assert p.fraction_corrected_ge1.between(-1e-9,1+1e-9).all() and p.fraction_raw_ge1.between(-1e-9,1+1e-9).all()
for age,g in e.groupby('age'):
 stage=p[p.grouping.eq('age')&p.age.eq(age)].set_index('gene');assert stage.nuclei.eq(len(g)).all()
 assert np.allclose(stage.loc[e.columns[3:],'mean_log1p_10k'],g.iloc[:,3:].mean(),atol=2e-5)
 for sample,lib in g.groupby('sample'):
  table=p[p.grouping.eq('age+sample')&p.age.eq(age)&p['sample'].eq(sample)].set_index('gene');assert table.nuclei.eq(len(lib)).all()
  assert np.allclose(table.loc[e.columns[3:],'mean_log1p_10k'],lib.iloc[:,3:].mean(),atol=2e-5)
a=pd.read_csv(OUT/'within_library_gene_associations.tsv',sep='\t');assert len(a)==48
for c in ['nucleus_spearman','QC_cycle_residual_spearman','raw_nucleus_spearman','cluster_mean_spearman']:assert a[c].between(-1,1).all()
ecm=pd.read_csv(OUT/'ECM_CM_coexpression_counts.tsv',sep='\t');assert ecm.nuclei.sum()==44463 and (ecm.ECM3_CM2<=ecm.ECM2_CM2).all() and (ecm.ECM2_CM2<=ecm.nuclei).all()
manifest={'status':'validated targeted literature-guided checks; no new subtype assignments','validation':['44,463 existing ventricular-enriched nuclei preserved in original order','One feature per selected gene; finite nonnegative expression','Stage/library profile means match per-cell data within float32/export tolerance','Detection bounds, association ranges and ECM accounting checked'],'scripts':[],'files':[]}
for name in ['check_ventricular_additional_axes.py','plot_ventricular_electrical_axis.py','ventricular_electrical_axis_sensitivity.py','validate_ventricular_additional_axes.py']:
 q=ROOT/'analysis'/name;manifest['scripts'].append({'file':str(q.relative_to(ROOT)),'sha256':hashlib.sha256(q.read_bytes()).hexdigest()})
for q in sorted(OUT.rglob('*')):
 if q.is_file() and q.name!='manifest.json':manifest['files'].append({'file':str(q.relative_to(OUT)),'bytes':q.stat().st_size,'sha256':hashlib.sha256(q.read_bytes()).hexdigest()})
(OUT/'manifest.json').write_text(json.dumps(manifest,indent=2));print(f'PASS:44,463 nuclei, {len(map)} genes, stage/library profiles, expression bounds and evidence checksums.')
