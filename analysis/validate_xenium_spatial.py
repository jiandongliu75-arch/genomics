"""Check raw-count/coordinate conservation and candidate-transfer rule invariants."""
from pathlib import Path
import hashlib
import json

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.neighbors import NearestNeighbors

ROOT = Path(__file__).resolve().parents[1]
BASE = Path('/workspace/heart-atlas')
OUT = ROOT / 'results/xenium_ventricular_mapping'
STAGES = ['P0', 'P7', 'P14', 'P21']
catalog = pd.read_csv(OUT / 'reference_group_catalog.tsv', sep='\t').set_index('microenvironment_cluster')
metadata = pd.read_csv(OUT / 'reference_metadata.tsv.gz', sep='\t', index_col=0)
reference = np.load(OUT / 'reference_shared_panel_counts.npz')
assert reference['corrected'].shape == reference['original'].shape == (73052, 381)
assert np.array_equal(reference['cell_ids'], metadata.index)
assert np.isfinite(reference['corrected']).all() and (reference['corrected'] >= 0).all()
assert (reference['original'] >= 0).all() and np.equal(reference['original'], np.floor(reference['original'])).all()
assert len(set(reference['genes'])) == 381 and 'Kcnd2' not in reference['genes'] and 'Irx5' not in reference['genes']
rng = np.random.default_rng(42)
samples, summaries = [], []


def vote(ref, query, labels):
    labels = np.asarray(labels, dtype=str)
    classes, codes = np.unique(labels, return_inverse=True)
    index = NearestNeighbors(n_neighbors=30, metric='cosine', algorithm='brute', n_jobs=4).fit(ref)
    distance, neighbor = index.kneighbors(query)
    weight = 1 / np.maximum(distance, .01)
    values = np.zeros((len(query), len(classes)), dtype=np.float32)
    for j in range(30):
        np.add.at(values, (np.arange(len(query)), codes[neighbor[:, j]]), weight[:, j])
    values /= values.sum(axis=1)[:, None]
    order = np.argsort(values, axis=1)
    return classes[order[:, -1]], values[np.arange(len(query)), order[:, -1]], values[np.arange(len(query)), order[:, -1]] - values[np.arange(len(query)), order[:, -2]]


for stage in STAGES:
    folder = OUT / stage
    data = ad.read_h5ad(folder / f'{stage}_xenium_counts_and_candidate_mapping.h5ad')
    obs = data.obs
    source = BASE / 'data/xenium' / str(obs.section.iloc[0])
    with h5py.File(source / 'cell_feature_matrix.h5') as handle:
        m = handle['matrix']
        original = sparse.csc_matrix((m['data'][:], m['indices'][:], m['indptr'][:]), shape=tuple(m['shape'][:])).T.tocsr()
        barcodes = np.array([v.decode() for v in m['barcodes'][:]])
        names = np.array([v.decode() for v in m['features/name'][:]])
        types = np.array([v.decode() for v in m['features/feature_type'][:]])
    assert data.n_vars == 403 and data.obs_names.is_unique and data.var_names.is_unique
    assert np.array_equal(barcodes, data.obs_names) and np.array_equal(names[types == 'Gene Expression'], data.var_names)
    assert (original[:, types == 'Gene Expression'] - data.X).nnz == 0
    assert np.issubdtype(data.X.dtype, np.integer) and (data.X.data >= 0).all()
    coords = pd.read_csv(source / 'cells.csv.gz').set_index('cell_id').loc[data.obs_names]
    assert np.allclose(obs[['x_centroid', 'y_centroid']], coords[['x_centroid', 'y_centroid']], atol=1e-9, rtol=0)
    assert np.array_equal(data.obsm['spatial'], coords[['x_centroid', 'y_centroid']].to_numpy())
    assert np.array_equal(np.asarray(data.X.sum(axis=1)).ravel(), obs.transcript_counts)
    for gene in ['Kcnd2', 'Irx5']:
        assert gene not in data.var_names and f'measured_{gene}_counts' not in obs
    assert 'NOT ASSAYED' in data.uns['Kcnd2_status']
    expected_qc = obs.transcript_counts.ge(30) & obs.measured_gene_count.ge(10) & obs.cell_area.gt(0) & obs.valid_coordinates & obs.negative_probe_fraction.le(.1)
    assert np.array_equal(expected_qc, obs.pass_expression_QC)
    assert np.array_equal(obs.missing_nucleus_area_flag, coords.nucleus_area.isna())
    overlap = obs.endothelial_regulators_detected.ge(1) & obs.endothelial_vascular_markers_detected.ge(3) & obs.positive_CM_fraction.ge(.10)
    assert np.array_equal(overlap, obs.CM_EC_program_overlap_review_flag)
    assert not (overlap & obs.measured_CM_marker_gate).any()
    myl2 = data.X[:, data.var_names.get_loc('Myl2')].toarray().ravel()
    myl4 = data.X[:, data.var_names.get_loc('Myl4')].toarray().ravel()
    expected_vent = obs.measured_CM_marker_gate.to_numpy() & (myl2 >= 2) & (myl2 >= myl4)
    vent = obs.spatial_lineage_candidate.eq('Ventricular-enriched CM').to_numpy()
    assert np.array_equal(expected_vent, vent)
    family = vent & obs.spatial_mapped_family.ne('Unresolved ventricular state').to_numpy()
    fine = obs.spatial_mapped_fine_group.ne('Not resolved').to_numpy()
    assert (~fine | family).all()
    assert obs.loc[family, 'spatial_family_score'].ge(.60 - 1e-6).all()
    assert obs.loc[family, 'shared_panel_transcripts'].ge(50).all()
    assert obs.loc[family, 'spatial_family_raw_reference_agreement'].all() and obs.loc[family, 'spatial_family_library_agreement'].all()
    assert obs.loc[fine, 'spatial_fine_group_score'].ge(.70 - 1e-6).all() and obs.loc[fine, 'shared_panel_transcripts'].ge(100).all()
    coarsecheck = pd.read_csv(folder / 'assay_aligned_family_recoverability.tsv', sep='\t').set_index('reference_class')
    finecheck = pd.read_csv(folder / 'assay_aligned_fine_group_recoverability.tsv', sep='\t').set_index('reference_class')
    assert coarsecheck.loc[obs.loc[family, 'spatial_mapped_family'], 'panel_supported'].all()
    assert finecheck.loc[obs.loc[fine, 'spatial_mapped_fine_group'], 'panel_supported'].all()
    assert catalog.loc[obs.loc[fine, 'spatial_mapped_fine_group'], 'both_library_raw_corrected_marker_count'].ge(3).all()
    assert np.array_equal(catalog.loc[obs.loc[fine, 'spatial_mapped_fine_group'], 'mapping_family'].to_numpy(), obs.loc[fine, 'spatial_mapped_family'].to_numpy())
    for field in ['X_assay_aligned_ventricular', 'X_assay_aligned_ventricular_original']:
        assert data.obsm[field].shape == (len(obs), 20)
        assert np.isfinite(data.obsm[field][vent]).all()
        assert np.isnan(data.obsm[field][~vent]).all()
    # Stratified checks reproduce votes from published reference/query embeddings,
    # independently of the code that assigned the labels.
    embedding = np.load(folder / 'reference_aligned_embeddings.npz')
    refobs = metadata.loc[embedding['ventricular_reference_cell_ids']]
    chosen = []
    for _, group in obs.loc[fine].groupby('spatial_mapped_fine_group', observed=True):
        ids = obs.index.get_indexer(group.index)
        chosen.extend(rng.choice(ids, min(30, len(ids)), replace=False))
    chosen = np.asarray(chosen, dtype=int)
    if len(chosen):
        target = data.obsm['X_assay_aligned_ventricular'][chosen]
        expected, score, margin = vote(embedding['ventricular_corrected'], target, refobs.ventricular_group)
        raw, _, _ = vote(embedding['ventricular_original'], data.obsm['X_assay_aligned_ventricular_original'][chosen], refobs.ventricular_group)
        library_predictions = []
        for library in sorted(refobs['sample'].unique()):
            selected = refobs['sample'].eq(library).to_numpy()
            prediction, _, _ = vote(embedding['ventricular_corrected'][selected], target, refobs.loc[selected, 'ventricular_group'])
            library_predictions.append(prediction)
        stored = obs.iloc[chosen].spatial_mapped_fine_group.to_numpy(dtype=str)
        assert np.array_equal(expected, stored) and np.array_equal(raw, stored)
        assert all(np.array_equal(value, stored) for value in library_predictions)
        assert (margin >= .20 - 1e-5).all() and np.allclose(score, obs.iloc[chosen].spatial_fine_group_score, atol=2e-5)
        for i, cell in enumerate(obs.index[chosen]):
            samples.append({'stage': stage, 'cell_id': cell, 'stored_group': stored[i], 'recomputed_corrected': expected[i], 'recomputed_original': raw[i], 'library1': library_predictions[0][i], 'library2': library_predictions[1][i], 'vote': float(score[i]), 'margin': float(margin[i])})
    image = json.loads((folder / 'morphology_display.json').read_text())
    assert image['channel'] == 'DAPI only'
    assert obs.x_centroid.between(0, image['extent_um'][1]).all() and obs.y_centroid.between(0, image['extent_um'][2]).all()
    for name in [f'{stage}_spatial_candidate_maps.png', f'{stage}_measured_marker_maps.png', f'{stage}_boundary_zoom.png']:
        assert (folder / name).stat().st_size > 10000
    summary = json.loads((folder / 'spatial_summary.json').read_text())
    assert summary['source_segmented_cells'] == len(obs) and summary['ventricular_CM_candidates'] == vent.sum() and summary['family_assigned'] == family.sum() and summary['fine_group_assigned'] == fine.sum()
    summaries.append(summary)
    print('VALIDATED_SECTION', stage, len(obs), 'unchanged count/coordinate rows;', len(chosen), 'fine-label vote checks', flush=True)

assert sum(row['source_segmented_cells'] for row in summaries) == 408054
pd.DataFrame(samples).to_csv(OUT / 'validated_fine_vote_samples.tsv', sep='\t', index=False)
manifest = {'status': 'Counts/coordinates and candidate-rule invariants validated; spatial biological accuracy not independently established.', 'validation': ['408,054 original segmented cells preserved across four sections; 403 measured gene targets each', 'Every original gene count and centroid conserved exactly against Zenodo inputs', 'QC, missing-nucleus flags, EC/CM overlap review and ventricular marker gates checked', 'Family/fine subset, score, depth, reference-support and source-marker invariants checked for all assigned cells', 'Stratified accepted fine-label votes, margins, original-reference and individual-library agreement independently recomputed from published embeddings', 'Physical image bounds, figure existence and summary accounting checked', 'Kcnd2/Irx5 absent; no measured or imputed gene-expression columns added'], 'fine_vote_samples': len(samples), 'files': [], 'scripts': []}
if (OUT / 'landmark_validation.json').exists():
    manifest['landmark_matching_validation'] = json.loads((OUT / 'landmark_validation.json').read_text())
for path in sorted(OUT.rglob('*')):
    if path.is_file() and path.name != 'manifest.json':
        assert path.stat().st_size < 95 * 1024 * 1024, f'GitHub file too large: {path}'
        manifest['files'].append({'file': str(path.relative_to(OUT)), 'bytes': path.stat().st_size, 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()})
for name in ['acquire_xenium_context.py', 'prepare_xenium_reference.py', 'map_xenium_ventricular.py', 'prepare_xenium_morphology.py', 'align_xenium_assays.py', 'finalize_xenium_spatial.py', 'plot_xenium_ventricular.py', 'check_xenium_landmark_matching.py', 'validate_xenium_landmarks.py', 'validate_xenium_spatial.py']:
    path = ROOT / 'analysis' / name
    manifest['scripts'].append({'file': str(path.relative_to(ROOT)), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()})
(OUT / 'manifest.json').write_text(json.dumps(manifest, indent=2))
print(f'PASS: 408,054 cells; exact original counts/coordinates; candidate rules; {len(samples)} fine-vote checks; figures/checksums.')
