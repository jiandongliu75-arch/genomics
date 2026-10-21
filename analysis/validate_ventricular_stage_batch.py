"""Validate population, stage preservation, graph/metadata alignment and count conservation."""
from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd,anndata as ad
ROOT=Path(__file__).resolve().parents[1];out=ROOT/'results/ventricular_stage_batch'
q=ad.read_h5ad(out/'ventricular_embeddings_metadata_only.h5ad')
a=pd.read_csv(ROOT/'results/cell_identities_subtypes/Cardiomyocyte/cell_subcluster_annotations.csv.gz',index_col=0);expected=a[a.proposed_subtype.eq('Ventricular-like')]
assert q.n_obs==43249 and q.n_vars==0 and q.obs_names.is_unique and q.obs_names.equals(expected.index)
assert q.obs.proposed_subtype.astype(str).eq('Ventricular-like').all()
assert q.obs.subcluster_id.astype(str).equals(expected.subcluster_id.astype(str))
representations=['uncorrected','within_stage_harmony','qc_sex_sensitivity']
for name in representations:
 z=q.obsm['X_'+name];assert z.shape==(43249,30) and np.isfinite(z).all()
 assert q.obsm['X_umap_'+name].shape==(43249,2) and np.isfinite(q.obsm['X_umap_'+name]).all()
 graph=q.obsp[name+'_connectivities'];assert graph.shape==(43249,43249) and graph.nnz>0
 for stage in ['P0','P7','P14','P21']:
  m=q.obs.age.eq(stage).to_numpy();assert np.allclose(z[m].mean(axis=0),q.obsm['X_uncorrected'][m].mean(axis=0),atol=1e-3)
ann=pd.read_csv(out/'ventricular_cells_embeddings_and_clusters.csv.gz',index_col=0);assert ann.index.equals(q.obs_names)
assert ann.ventricular_cluster.equals(q.obs.ventricular_cluster.astype(str))
s=pd.read_csv(out/'ventricular_cluster_summary.tsv',sep='\t');assert s.nuclei.sum()==43249 and s.ventricular_cluster.is_unique
pb=pd.read_csv(out/'ventricular_cluster_corrected_pseudobulk.tsv.gz',sep='\t',index_col=0);assert len(pb)==32285 and set(pb.columns[1:])==set(q.obs.ventricular_cluster)
assert np.isclose(pb.iloc[:,1:].to_numpy().sum(),expected.total_counts.sum(),rtol=1e-5), 'Corrected counts must remain unchanged by embedding adjustment'
mix=pd.read_csv(out/'within_stage_sample_mixing.tsv',sep='\t');assert len(mix)==12
assert mix.different_library_neighbor_fraction.between(0,1).all()
assert mix.mean_normalized_library_entropy.between(-1e-8,1+1e-8).all()
manifest={'status':'validated','validation':['43,249 previously ventricular-like nuclei, unique and aligned to original annotations','All three PCA/UMAP representations and graphs have aligned, finite coordinates','Stage PCA centroids preserved in primary and sensitivity representations','Original corrected counts conserved in all-gene pseudobulk export','All cluster and per-stage metric counts accounted for'],'script_sha256':hashlib.sha256((ROOT/'analysis/ventricular_stage_batch.py').read_bytes()).hexdigest(),'scripts':[{'file':str(f.relative_to(ROOT)),'sha256':hashlib.sha256(f.read_bytes()).hexdigest()} for f in [ROOT/'analysis/ventricular_stage_batch.py',ROOT/'analysis/ventricular_within_stage_markers.py',ROOT/'analysis/validate_ventricular_stage_batch.py']],'files':[]}
for p in sorted(out.iterdir()):
 if p.is_file() and p.name!='manifest.json':manifest['files'].append({'name':p.name,'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
assert all(f['bytes']<90_000_000 for f in manifest['files'])
(out/'manifest.json').write_text(json.dumps(manifest,indent=2));print('PASS: ventricular population, coordinates/graphs, preserved stage centroids, original annotation alignment, count conservation, metrics and checksums.')
