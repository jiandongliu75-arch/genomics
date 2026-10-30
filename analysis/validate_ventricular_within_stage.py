"""Validate independently resolved ventricular stages, counts and support annotations."""
from pathlib import Path
import hashlib,json
import numpy as np,pandas as pd,anndata as ad
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/ventricular_within_stage'
ann=pd.read_csv(OUT/'cell_annotations.csv.gz',index_col=0);summary=pd.read_csv(OUT/'cluster_summary.tsv',sep='\t');source=pd.read_csv(ROOT/'results/ventricular_stage_batch/ventricular_cells_embeddings_and_clusters.csv.gz',index_col=0)
assert len(ann)==43249 and ann.index.is_unique and set(ann.index)==set(source.index)
assert ann.proposed_subtype.eq('Ventricular-like').all() and ann.stage_state_label.notna().all() and ann.stage_support_status.notna().all()
assert summary.stage_subcluster_id.is_unique and summary.nuclei.sum()==43249
for stage in ['P0','P7','P14','P21']:
 q=ad.read_h5ad(OUT/stage/'embeddings_metadata_only.h5ad');s=ann[ann.age.eq(stage)];assert q.n_vars==0 and q.n_obs==len(s) and q.obs_names.equals(s.index)
 assert q.obs.stage_subcluster_id.astype(str).equals(s.stage_subcluster_id.astype(str))
 assert q.obs.stage_state_label.astype(str).equals(s.stage_state_label.astype(str))
 assert np.isfinite(q.obsm['X_umap']).all() and np.allclose(q.obsm['X_umap'],s[['within_stage_UMAP1','within_stage_UMAP2']].to_numpy(),atol=1e-5)
 for key in ['uncorrected','harmony','qc_sex_sensitivity']:
  assert q.obsm['X_'+key].shape==(len(s),30) and np.isfinite(q.obsm['X_'+key]).all();assert q.obsp[key+'_connectivities'].shape==(len(s),len(s))
  assert np.allclose(q.obsm['X_'+key].mean(axis=0),q.obsm['X_uncorrected'].mean(axis=0),atol=1e-3)
 pb=pd.read_csv(OUT/stage/'corrected_pseudobulk.tsv.gz',sep='\t',index_col=0);assert len(pb)==32285 and set(pb.columns[1:])==set(s.stage_subcluster_id)
 assert np.isclose(pb.iloc[:,1:].to_numpy().sum(),source.loc[s.index,'total_counts'].sum(),rtol=1e-5),stage+' count conservation'
 for cid,g in s.groupby('stage_subcluster_id'):
  row=summary.set_index('stage_subcluster_id').loc[cid];assert row.nuclei==len(g) and g.stage_support_status.eq(row.stage_support_status).all()
rep=summary.replicate_program_support.astype(bool)
assert (summary.loc[rep,'nuclei']>=100).all() and (summary.loc[rep,'libraries_ge30']==2).all() and (summary.loc[rep,'largest_library_fraction']<=.9).all()
assert (summary.loc[rep,'both_library_marker_count']>=3).all() and (summary.loc[rep,'seed_min_best_jaccard']>=.5).all() and (summary.loc[rep,'seed_mean_best_jaccard']>=.65).all()
passed=summary.stage_support_status.eq('supported transcriptional state');assert (summary.loc[passed,'qc_sex_sensitive']==False).all() and (summary.loc[passed,'sample_adjustment_sensitive']==False).all()
markers=pd.read_csv(OUT/'within_stage_markers.tsv',sep='\t');assert (markers.positive_libraries<=markers.eligible_libraries).all()
assert markers.fraction_corrected_ge1_cluster.between(-1e-9,1+1e-9).all()
mix=pd.read_csv(OUT/'sample_mixing.tsv',sep='\t');assert len(mix)==12 and mix.different_library_neighbor_fraction.between(0,1).all()
manifest={'status':'validated independent within-stage analysis','validation':['43,249 original ventricular nuclei preserved, unique and aligned','Independent stage coordinates/graphs, CSV and metadata H5AD align','Original corrected count totals conserved in each stage pseudobulk','Support labels meet prespecified marker/library/seed/nuisance heuristics','Detection fractions, neighbor metrics and artifact sizes checked'],'scripts':[],'files':[]}
for name in ['ventricular_within_stage.py','finalize_ventricular_within_stage_labels.py','plot_ventricular_within_stage.py','validate_ventricular_within_stage.py']:
 p=ROOT/'analysis'/name;manifest['scripts'].append({'file':str(p.relative_to(ROOT)),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
for p in sorted(OUT.rglob('*')):
 if p.is_file() and p.name!='manifest.json':manifest['files'].append({'file':str(p.relative_to(OUT)),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
assert all(f['bytes']<90_000_000 for f in manifest['files'])
(OUT/'manifest.json').write_text(json.dumps(manifest,indent=2));print(f'PASS: 43,249 nuclei; {len(summary)} within-stage groups; preserved counts; aligned coordinates/graphs; support gates; checksums.')
