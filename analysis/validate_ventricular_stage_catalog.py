"""Check catalog accounting and reproduce similarity from exported signatures."""
from pathlib import Path
import hashlib
import json

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/ventricular_stage_subtype_catalog'
SOURCE = ROOT / 'results/ventricular_microenvironment'
catalog = pd.read_csv(OUT / 'stage_subtype_catalog.tsv', sep='\t')
old = pd.read_csv(SOURCE / 'cluster_summary.tsv', sep='\t')
assert len(catalog) == 28 and catalog.nuclei.sum() == 44463
for field in old.columns:
    pd.testing.assert_series_equal(catalog[field], old[field])
assert catalog.candidate_identity.notna().all() and catalog.potential_function_and_limits.notna().all()
assert np.allclose(catalog.groupby('stage').fraction_within_retained_stage_CMs.sum(), 1)
cells = pd.read_csv(OUT / 'cell_catalog_and_comparison_memberships.tsv.gz', sep='\t', index_col=0)
original = pd.read_csv(SOURCE / 'cell_annotations_and_scores.csv.gz', index_col=0)
assert len(cells) == 44463 and cells.index.equals(original.index)
for field in ['age', 'sample', 'microenvironment_cluster', 'gate_origin']:
    pd.testing.assert_series_equal(cells[field], original[field])
assert cells.groupby('microenvironment_cluster').size().sort_index().equals(catalog.set_index('microenvironment_cluster').nuclei.sort_index().rename(None))
lookup = catalog.set_index('microenvironment_cluster')
for field in ['candidate_identity', 'functional_family', 'support_status']:
    assert cells[field].equals(cells.microenvironment_cluster.map(lookup[field]).rename(field))
comparison = pd.read_csv(OUT / 'cross_stage_comparison_sets.tsv', sep='\t').fillna('')
for stage in ['P0', 'P7', 'P14', 'P21']:
    valid = set(catalog.loc[catalog.stage.eq(stage), 'microenvironment_cluster'])
    for value in comparison[stage]:
        assert set(filter(None, value.split(';'))) <= valid
features = pd.read_csv(OUT / 'signature_feature_mapping.tsv', sep='\t', index_col=0)
assert features.gene_symbol.is_unique and features.index.is_unique
assert (features.interpreted_gene_excluded_similarity_feature <= features.primary_similarity_feature).all()
signatures = pd.read_csv(OUT / 'within_stage_relative_aggregate_signatures.tsv.gz', sep='\t', index_col=0)
assert set(signatures.columns) == set(catalog.microenvironment_cluster)
assert signatures.index.equals(features.index[features.primary_similarity_feature])
assert np.isfinite(signatures.to_numpy()).all()
primary = pd.read_csv(OUT / 'relative_signature_similarity_matrix.tsv', sep='\t', index_col=0)
excluded = pd.read_csv(OUT / 'interpreted_genes_excluded_similarity_matrix.tsv', sep='\t', index_col=0)
recomputed = signatures.corr(method='spearman')
assert np.allclose(primary.loc[recomputed.index, recomputed.columns], recomputed, atol=2e-5)
recomputed_excluded = signatures.loc[features.index[features.interpreted_gene_excluded_similarity_feature]].corr(method='spearman')
assert np.allclose(excluded.loc[recomputed_excluded.index, recomputed_excluded.columns], recomputed_excluded, atol=2e-5)
for matrix in [primary, excluded]:
    assert np.allclose(matrix, matrix.T) and np.allclose(np.diag(matrix), 1)
    assert np.abs(matrix.to_numpy()).max() <= 1 + 1e-9
pairs = pd.read_csv(OUT / 'cross_stage_pairwise_similarity.tsv', sep='\t')
assert len(pairs) == 578 and not pairs.duplicated(['source_cluster', 'target_cluster']).any()
assert pairs.source_stage.ne(pairs.target_stage).all()
for row in pairs.itertuples():
    assert np.isclose(row.relative_signature_spearman, primary.loc[row.source_cluster, row.target_cluster])
    assert np.isclose(row.interpreted_genes_excluded_spearman, excluded.loc[row.source_cluster, row.target_cluster])
for (source, stage), group in pairs.groupby(['source_cluster', 'target_stage']):
    assert group.rank_within_target_stage.min() == 1
    ranks = group.relative_signature_spearman.rank(ascending=False, method='min').to_numpy()
    assert np.array_equal(ranks, group.rank_within_target_stage)
manifest = {'validation': ['28 original groups and 44,463 nuclei preserved with all prior support fields', 'Candidate functions and overlapping proposed comparison sets refer only to valid same-stage groups', 'Per-cell annotations preserve source membership and map to catalog', 'One feature per symbol; finite signatures and nested exclusion feature set', 'Primary and interpreted-gene-excluded similarity matrices reproduce from exported signatures', 'All 578 directed cross-stage comparisons align with matrices and within-stage ranks'], 'files': [], 'scripts': []}
for path in sorted(OUT.rglob('*')):
    if path.is_file() and path.name != 'manifest.json':
        manifest['files'].append({'file': str(path.relative_to(OUT)), 'bytes': path.stat().st_size, 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()})
for name in ['ventricular_stage_subtype_catalog.py', 'plot_ventricular_stage_comparisons.py', 'validate_ventricular_stage_catalog.py']:
    path = ROOT / 'analysis' / name
    manifest['scripts'].append({'file': str(path.relative_to(ROOT)), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()})
(OUT / 'manifest.json').write_text(json.dumps(manifest, indent=2))
print('PASS: 28 original groups / 44,463 nuclei, per-cell catalog, comparison membership and all 578 cross-stage similarities.')
