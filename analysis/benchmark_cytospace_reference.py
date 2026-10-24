"""Technical shared-panel recoverability across the two snRNA libraries per age.

The queries are nuclei, not Xenium cells: this is not spatial accuracy validation.
Stage groups were defined using both libraries, so these are technical checks of
panel/algorithm recoverability rather than independent biological labels.
"""
from pathlib import Path
import gzip
import json
import os
import shutil
import subprocess

import numpy as np
import pandas as pd
from sklearn.metrics import precision_recall_fscore_support

from run_cytospace_xenium import SOURCE, OUT, CACHE, CYTO, BASE, STAGES, REF_LINEAGES, write_matrix


def main():
    meta = pd.read_csv(SOURCE / 'reference_metadata.tsv.gz', sep='\t', index_col=0)
    panel = np.load(SOURCE / 'reference_shared_panel_counts.npz')
    counts = panel['corrected']
    valid = meta.broad_class.isin(REF_LINEAGES) & (counts.sum(1) >= 20) & ((counts > 0).sum(1) >= 8)
    # Use the same reference eligibility as the spatial runs.
    raw = panel['original']
    valid &= (raw.sum(1) >= 20) & ((raw > 0).sum(1) >= 8)
    all_metrics = []
    for stage in STAGES:
        libraries = sorted(meta.loc[valid & meta.age.eq(stage), 'sample'].unique())
        for fold, query_library in enumerate(libraries):
            dest = OUT / 'reference_panel_benchmark' / stage / f'query_library{fold+1}'
            dest.mkdir(parents=True, exist_ok=True)
            work = CACHE / 'benchmark' / stage / f'library{fold+1}'
            qi = np.flatnonzero(valid & meta.age.eq(stage) & meta['sample'].eq(query_library))
            ri = np.flatnonzero(valid & meta.age.eq(stage) & ~meta['sample'].eq(query_library))
            if not (dest / 'assigned_locations.csv.gz').exists():
                rng = np.random.default_rng(42)
                qi, ri = rng.permutation(qi), rng.permutation(ri)
                qtypes = meta.iloc[qi].broad_class.map(REF_LINEAGES)
                rtypes = meta.iloc[ri].broad_class.map(REF_LINEAGES)
                common = set(qtypes) & set(rtypes)
                qi, ri = qi[qtypes.isin(common)], ri[rtypes.isin(common)]
                write_matrix(work / 'reference', counts[ri], panel['genes'], meta.index[ri])
                write_matrix(work / 'query', counts[qi], panel['genes'], meta.index[qi])
                meta.iloc[ri].broad_class.map(REF_LINEAGES).to_frame('CellType').to_csv(work / 'reference_labels.tsv', sep='\t', index_label='CellID')
                meta.iloc[qi].broad_class.map(REF_LINEAGES).to_frame('CellType').to_csv(work / 'query_labels.tsv', sep='\t', index_label='SpotID')
                pd.DataFrame({'x_centroid': np.arange(len(qi)), 'y_centroid': np.zeros(len(qi))}, index=meta.index[qi]).to_csv(work / 'coordinates.tsv', sep='\t', index_label='SpotID')
                cmd = [str(CYTO), '--single-cell', '--scRNA-path', str(work / 'reference/matrix.mtx.gz'),
                       '--cell-type-path', str(work / 'reference_labels.tsv'), '--st-path', str(work / 'query/matrix.mtx.gz'),
                       '--coordinates-path', str(work / 'coordinates.tsv'), '--st-cell-type-path', str(work / 'query_labels.tsv'),
                       '--output-folder', str(work / 'upstream_output'), '--number-of-selected-spots', '3000',
                       '--number-of-processors', '2', '--solver-method', 'lapjv', '--seed', '42', '--downsample-off', '--plot-off']
                env = {**os.environ, 'OPENBLAS_NUM_THREADS': '1', 'OMP_NUM_THREADS': '1', 'MKL_NUM_THREADS': '1',
                       'NUMBA_CACHE_DIR': '/tmp/cytospace-numba', 'MPLCONFIGDIR': '/tmp/heart-mpl', 'PYTHONUNBUFFERED': '1'}
                print(f'BENCHMARK {stage} query={query_library} n={len(qi)}', flush=True)
                with open(dest / 'upstream_stdout.txt', 'w') as log:
                    subprocess.run(cmd, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
                assigned = pd.read_csv(work / 'upstream_output/assigned_locations.csv')
                assigned.to_csv(dest / 'assigned_locations.csv.gz', index=False)
                shutil.copyfile(work / 'upstream_output/log.txt', dest / 'upstream_log.txt')
                (dest / 'benchmark_manifest.json').write_text(json.dumps({'stage': stage, 'query_library': query_library,
                    'training_library': [x for x in libraries if x != query_library], 'command': cmd, 'genes': len(panel['genes']),
                    'query_counts': 'DecontX-corrected count estimates; no Xenium technology effects simulated',
                    'coordinates': 'arbitrary index and zero; no tissue locations, never used in assignment cost',
                    'limitation': 'Technical cross-library panel recoverability of existing clusters; not independent Xenium accuracy'}, indent=2))
                shutil.rmtree(work)
            assigned = pd.read_csv(dest / 'assigned_locations.csv.gz')
            assert assigned.SpotID.is_unique
            truth = meta.loc[assigned.SpotID]
            predicted = meta.loc[assigned.OriginalCID]
            assert truth['sample'].eq(query_library).all()
            assert not predicted['sample'].eq(query_library).any()
            vent = truth.broad_class.eq('Ventricular-enriched CM').to_numpy()
            for field in ['ventricular_group', 'mapping_family']:
                y = truth.iloc[np.flatnonzero(vent)][field].astype(str).to_numpy()
                yh = predicted.iloc[np.flatnonzero(vent)][field].astype(str).to_numpy()
                classes = sorted(set(y) | set(yh))
                precision, recall, f1, n = precision_recall_fscore_support(y, yh, labels=classes, zero_division=0)
                for label, p, r, f, qn in zip(classes, precision, recall, f1, n):
                    tn = (meta.iloc[ri][field].eq(label)).sum()
                    all_metrics.append({'stage': stage, 'query_library': query_library, 'resolution': field,
                        'label': label, 'precision': p, 'recall': r, 'f1': f, 'query_nuclei': int(qn), 'training_nuclei': int(tn)})
                pd.crosstab(pd.Series(y, name='known_reference_group'), pd.Series(yh, name='assigned_reference_group')).to_csv(dest / f'{field}_confusion.tsv', sep='\t')
    metrics = pd.DataFrame(all_metrics)
    metrics.to_csv(OUT / 'reference_panel_benchmark/heldout_library_metrics.tsv', sep='\t', index=False)
    summary = metrics.groupby(['stage', 'resolution', 'label']).agg(min_precision=('precision', 'min'),
        min_recall=('recall', 'min'), min_f1=('f1', 'min'), min_query_nuclei=('query_nuclei', 'min'),
        min_training_nuclei=('training_nuclei', 'min'), tested_folds=('f1', 'size')).reset_index()
    summary['panel_recoverable'] = (summary.tested_folds.eq(2) & (summary.min_precision >= .60) &
        (summary.min_recall >= .35) & (summary.min_f1 >= .50) & (summary.min_query_nuclei >= 30) & (summary.min_training_nuclei >= 30))
    summary.to_csv(OUT / 'reference_panel_benchmark/panel_recoverability.tsv', sep='\t', index=False)
    sensitivity = []
    for cutoff in [.50, .55, .60, .65, .70]:
        table = summary.copy()
        table['minimum_precision_cutoff'] = cutoff
        table['passes_at_this_cutoff'] = (table.tested_folds.eq(2) & (table.min_precision >= cutoff) &
            (table.min_recall >= .35) & (table.min_f1 >= .50) & (table.min_query_nuclei >= 30) & (table.min_training_nuclei >= 30))
        sensitivity.append(table)
    pd.concat(sensitivity, ignore_index=True).to_csv(OUT / 'reference_panel_benchmark/precision_threshold_sensitivity.tsv', sep='\t', index=False)
    print(summary.loc[summary.panel_recoverable].to_string(index=False), flush=True)


if __name__ == '__main__':
    main()
