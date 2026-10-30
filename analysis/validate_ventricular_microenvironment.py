"""Population, score, gene-feature, count-preservation and evidence validation."""
from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd,anndata as ad
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/ventricular_microenvironment'
cfg=json.loads((OUT/'program_definitions.json').read_text());programs=cfg['programs'];ann=pd.read_csv(OUT/'cell_annotations_and_scores.csv.gz',index_col=0);s=pd.read_csv(OUT/'cluster_summary.tsv',sep='\t');base=pd.read_csv(ROOT/'results/cell_identities_subtypes/all_cell_identity_and_subtype_annotations.csv.gz',index_col=0)
expected=base[base.proposed_subtype.eq('Ventricular-like')|base.subcluster_id.isin(['CAR_3','CAR_9'])]
assert len(ann)==44463 and ann.index.is_unique and set(ann.index)==set(expected.index)
assert ann.gate_origin.eq('core_ventricular').sum()==43249 and ann.gate_origin.eq('Slit2_overlap_reference').sum()==1214
assert set(ann.loc[ann.gate_origin.eq('Slit2_overlap_reference'),'subcluster_id'])=={'CAR_3','CAR_9'}
assert ann.microenvironment_label.notna().all() and ann.support_status.notna().all() and s.nuclei.sum()==44463 and s.microenvironment_cluster.is_unique
mapping=pd.read_csv(OUT/'program_gene_mapping.tsv',sep='\t').fillna('');assert mapping.present.all()
assert mapping.score_gene_ids.map(lambda x:len(x.split(';'))==1).all(), 'One scoring feature per requested symbol'
for stage in ['P0','P7','P14','P21']:
 q=ad.read_h5ad(OUT/stage/'embeddings_and_scores_metadata_only.h5ad');a=ann[ann.age.eq(stage)]
 assert q.n_vars==0 and q.n_obs==len(a) and q.obs_names.equals(a.index)
 assert q.obs.microenvironment_cluster.astype(str).equals(a.microenvironment_cluster.astype(str))
 assert np.allclose(q.obsm['X_umap'],a[['UMAP1','UMAP2']].to_numpy(),atol=1e-5) and np.isfinite(q.obsm['X_umap']).all()
 for name in ['all_gene_HVG','uncorrected','hypothesis_genes_excluded','QC_sex_cycle_sensitivity']:
  assert q.obsm['X_'+name].shape==(len(a),30) and np.isfinite(q.obsm['X_'+name]).all()
  assert q.obsp[name+'_connectivities'].shape==(len(a),len(a))
 for key in programs:
  for prefix in ['score_','z_','raw_z_','qc_adjusted_z_','coverage_','raw_coverage_']:assert np.isfinite(a[prefix+key]).all()
  for prefix in ['coverage_','raw_coverage_']:assert a[prefix+key].between(-1e-9,1+1e-9).all()
  for lib,g in a.groupby('sample'):
   for prefix in ['z_','raw_z_','qc_adjusted_z_']:
    vals=g[prefix+key].to_numpy();assert abs(vals.mean())<2e-5 and (np.isclose(vals.std(),1,atol=2e-5) or vals.std()<1e-8)
 pb=pd.read_csv(OUT/stage/'corrected_pseudobulk.tsv.gz',sep='\t',index_col=0);assert len(pb)==32285 and set(pb.columns[1:])==set(a.microenvironment_cluster)
 assert np.isclose(pb.iloc[:,1:].to_numpy().sum(),expected.loc[a.index,'total_counts'].sum(),rtol=1e-5),stage+' original corrected-count conservation'
 for cid,g in a.groupby('microenvironment_cluster'):
  row=s.set_index('microenvironment_cluster').loc[cid];assert row.nuclei==len(g)
  assert g.support_status.eq(row.support_status).all() and g.microenvironment_label.eq(row.microenvironment_label).all()
rep=s.support_status.ne('candidate/review');assert (s.loc[rep,'nuclei']>=100).all() and (s.loc[rep,'libraries_ge30']==2).all() and (s.loc[rep,'largest_library_fraction']<=.9).all()
assert (s.loc[rep,'both_library_raw_corrected_marker_count']>=3).all() and (s.loc[rep,'seed_mean_best_jaccard']>=.65).all() and (s.loc[rep,'seed_min_best_jaccard']>=.5).all()
supported=s.support_status.eq('supported RNA state')
for col,minimum in [('QC_sensitivity_best_jaccard',.5),('panel_excluded_best_jaccard',.5),('uncorrected_best_jaccard',.4)]:assert (s.loc[supported,col]>=minimum).all()
e=pd.read_csv(OUT/'cluster_program_evidence.tsv',sep='\t').fillna('');passed=e[e.program_rule_pass]
for _,r in passed.iterrows():
 definition=programs[r.program];symbols=r.supported_genes.split(';') if r.supported_genes else []
 assert len(set(symbols))==r.supported_gene_count and r.supported_gene_count>=definition['min_supported_genes']
 assert r.corrected_z_effect>.35 and r.raw_z_effect>0
 assert min(map(float,r.library_corrected_effects.split(';')))>.15 and min(map(float,r.library_raw_effects.split(';')))>0
for _,r in s[s.selected_program.eq('HIF_response')].iterrows():assert bool(e[(e.cluster==r.microenvironment_cluster)&e.program.eq('HIF_non_glycolytic')].program_rule_pass.iloc[0])
quality=pd.read_csv(OUT/'module_score_quality.tsv',sep='\t');assert len(quality)==4*len(programs)
cohort=pd.read_csv(OUT/'cohort_summary.tsv',sep='\t');assert cohort.expanded_nuclei.sum()==44463 and cohort.core_nuclei.sum()==43249 and cohort.reference_nuclei.sum()==1214
manifest={'status':'validated','validation':['44,463 unique expanded-gate nuclei =43,249 core +1,214 tracked ambiguous references','Independent stage UMAP/PC/graph/CSV alignment','Per-library score centering/scaling, finite values and detection coverage','One scoring feature per requested gene symbol; distinct-symbol support','Original corrected counts conserved; metadata-only files have no expression','State/program support gates and HIF specificity guard','Source provenance and artifact hashes'],'scripts':[],'files':[]}
for name in ['ventricular_microenvironment.py','plot_ventricular_microenvironment.py','validate_ventricular_microenvironment.py','finalize_ventricular_microenvironment.py']:
 p=ROOT/'analysis'/name
 if p.exists():manifest['scripts'].append({'file':str(p.relative_to(ROOT)),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
for p in sorted(OUT.rglob('*')):
 if p.is_file() and p.name!='manifest.json':manifest['files'].append({'file':str(p.relative_to(OUT)),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
assert all(f['bytes']<90_000_000 for f in manifest['files'])
(OUT/'manifest.json').write_text(json.dumps(manifest,indent=2));print(f'PASS:44,463 nuclei; {len(s)} clusters; core/reference gate; scores/features; preserved counts; state/program gates; checksums.')
