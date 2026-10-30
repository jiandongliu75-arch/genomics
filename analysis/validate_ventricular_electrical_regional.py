"""Validate published targeted matrices and within-library contrasts."""
from pathlib import Path
import hashlib
import json

import numpy as np
import pandas as pd
from scipy import sparse

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/ventricular_electrical_regional'
ann = pd.read_csv(OUT / 'cell_metadata.tsv.gz', sep='\t', index_col=0)
source = pd.read_csv(ROOT / 'results/ventricular_microenvironment/cell_annotations_and_scores.csv.gz', index_col=0)
assert len(ann) == 44463 and ann.index.is_unique and ann.index.equals(source.index)
assert ann.equals(source[ann.columns])
mapping = pd.read_csv(OUT / 'gene_mapping.tsv', sep='\t')
assert len(mapping) == 17 and mapping.gene.is_unique and mapping.gene_id.is_unique
totals = pd.read_csv(OUT / 'normalization_totals.tsv.gz', sep='\t', index_col=0)
expr = pd.read_csv(OUT / 'targeted_log_expression.csv.gz', index_col=0)
assert expr.index.equals(ann.index) and totals.index.equals(ann.index)
assert list(expr.columns) == list(mapping.gene)
original = sparse.load_npz(OUT / 'targeted_original_counts.npz').toarray()
corrected = sparse.load_npz(OUT / 'targeted_corrected_counts.npz').toarray()
assert original.shape == corrected.shape == (44463, 17)
assert (original >= 0).all() and np.equal(original, np.floor(original)).all()
assert np.isfinite(corrected).all() and (corrected >= 0).all()
assert (totals > 0).all().all()
logs = np.log1p(corrected * 10000 / totals.corrected_whole_nucleus_total.to_numpy()[:, None])
rawlogs = np.log1p(original * 10000 / totals.original_whole_nucleus_total.to_numpy()[:, None])
assert np.allclose(expr, logs, atol=2e-5)
profiles = pd.read_csv(OUT / 'gene_profiles.tsv', sep='\t')
for age, idx in ann.groupby('age').indices.items():
    prof = profiles[profiles.grouping.eq('age') & profiles.age.eq(age)].set_index('gene').loc[mapping.gene]
    assert prof.nuclei.eq(len(idx)).all()
    assert np.allclose(prof.mean_log1p_10k, logs[idx].mean(axis=0), atol=2e-5)
    assert np.allclose(prof.mean_raw_log1p_10k, rawlogs[idx].mean(axis=0), atol=2e-5)
    assert np.allclose(prof.fraction_corrected_ge1, (corrected[idx] >= 1).mean(axis=0))
    assert np.allclose(prof.fraction_raw_ge1, (original[idx] >= 1).mean(axis=0))
contrasts = pd.read_csv(OUT / 'cluster_gene_contrasts.tsv', sep='\t')
assert len(contrasts) == 28 * 17 and not contrasts.duplicated(['cluster', 'gene']).any()
assert contrasts.groupby('cluster').nuclei.first().sum() == len(ann)
for (age, cluster), sub in contrasts.groupby(['age', 'cluster']):
    expected, rawexpected = [], []
    for sample in ann.loc[ann.age.eq(age), 'sample'].unique():
        lib = ann.age.eq(age).to_numpy() & ann['sample'].eq(sample).to_numpy()
        inside = np.flatnonzero(lib & ann.microenvironment_cluster.eq(cluster).to_numpy())
        outside = np.flatnonzero(lib & ~ann.microenvironment_cluster.eq(cluster).to_numpy())
        if len(inside) >= 30 and len(outside) >= 30:
            expected.append(logs[inside].mean(axis=0) - logs[outside].mean(axis=0))
            rawexpected.append(rawlogs[inside].mean(axis=0) - rawlogs[outside].mean(axis=0))
    tab = sub.set_index('gene').loc[mapping.gene]
    assert tab.eligible_libraries.eq(len(expected)).all()
    assert np.allclose(tab.corrected_delta, np.mean(expected, axis=0), atol=2e-5)
    assert np.allclose(tab.raw_delta, np.mean(rawexpected, axis=0), atol=2e-5)
for name in ['ventricular_electrical_regional_contrasts.png', 'ventricular_electrical_regional_contrasts.pdf']:
    assert (OUT / name).stat().st_size > 10000
manifest = {'validation': ['44,463 original nuclei and annotations unchanged', '17-feature count matrices, integer original counts, nonnegative finite corrected counts', 'Whole-nucleus normalization reconstructs exported expression', 'Stage expression/detection profiles reconstruct from original and corrected counts', 'All 476 expanded cluster/gene contrasts reconstruct independently from published matrices and library membership'], 'files': [], 'scripts': []}
for path in sorted(OUT.rglob('*')):
    if path.is_file() and path.name != 'manifest.json':
        manifest['files'].append({'file': str(path.relative_to(OUT)), 'bytes': path.stat().st_size, 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()})
for name in ['ventricular_electrical_regional.py', 'plot_ventricular_electrical_regional.py', 'validate_ventricular_electrical_regional.py']:
    path = ROOT / 'analysis' / name
    manifest['scripts'].append({'file': str(path.relative_to(ROOT)), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()})
(OUT / 'manifest.json').write_text(json.dumps(manifest, indent=2))
print('PASS: 44,463 unchanged nuclei, 17-gene original/corrected matrices, normalization, stage profiles and all 476 expanded contrasts.')
