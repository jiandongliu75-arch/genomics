"""Default CytoSPACE parameter comparison on the published ventricular cohort.

Published final spatial states are a comparison layer, not supplied subtypes.
The snRNA reference is independently QC-reviewed/annotated, so this is not an
exact reproduction of the authors' unavailable input objects/label files.
"""
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import time

import anndata as ad
import numpy as np
import pandas as pd

from run_cytospace_xenium import SOURCE, BASE, CYTO, write_matrix, lognorm, paired_r, digest
from acquire_author_heart_atlas import OUT, STAGES

CACHE = BASE / 'results/author_parameter_baseline_cache'
VENT_STATES = ['v.CM', 'Ankrd1.CM', 'Slit2.CM', 'Ddc.CM', 'p.CM']


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--stages', nargs='+', default=STAGES, choices=STAGES)
    args = parser.parse_args()
    metadata = pd.read_csv(SOURCE / 'reference_metadata.tsv.gz', sep='\t', index_col=0)
    bundle = np.load(SOURCE / 'reference_shared_panel_counts.npz')
    raw = bundle['original']
    assert np.array_equal(raw, raw.astype(np.int32))
    raw = raw.astype(np.int32)
    genes = bundle['genes'].astype(str)
    for stage in args.stages:
        folder = OUT / stage / 'independent_default_parameters'
        if (folder / 'run_manifest.json').exists():
            print('SKIP_COMPLETED', stage, flush=True)
            continue
        folder.mkdir(parents=True, exist_ok=True)
        labels = pd.read_csv(OUT / stage / 'original_cells_with_author_labels.tsv.gz', sep='\t', index_col=0, float_precision='round_trip')
        data = ad.read_h5ad(SOURCE / stage / f'{stage}_xenium_counts_and_candidate_mapping.h5ad')
        assert np.array_equal(labels.index, data.obs_names)
        query = data[:, genes].X.toarray().astype(np.int32)
        included = labels.author_atlas_cell_state.isin(VENT_STATES) & (query.sum(1) > 0)
        ri = np.flatnonzero(metadata.age.eq(stage) & metadata.broad_class.eq('Ventricular-enriched CM') & (raw.sum(1) > 0))
        qi = np.flatnonzero(included)
        rng = np.random.default_rng(1)
        ri, qi = rng.permutation(ri), rng.permutation(qi)
        reference = metadata.iloc[ri]
        work = CACHE / stage
        write_matrix(work / 'reference', raw[ri], genes, reference.index)
        write_matrix(work / 'query', query[qi], genes, labels.index[qi])
        pd.DataFrame({'CellType': 'Ventricular_CM'}, index=reference.index).to_csv(work / 'reference_labels.tsv', sep='\t', index_label='CellID')
        pd.DataFrame({'CellType': 'Ventricular_CM'}, index=labels.index[qi]).to_csv(work / 'query_labels.tsv', sep='\t', index_label='SpotID')
        labels.iloc[qi][['x_centroid', 'y_centroid']].to_csv(work / 'coordinates.tsv', sep='\t', index_label='SpotID')
        reference.to_csv(folder / 'reference_pool.tsv.gz', sep='\t')
        command = [str(CYTO), '--single-cell', '--scRNA-path', str(work / 'reference/matrix.mtx.gz'),
            '--cell-type-path', str(work / 'reference_labels.tsv'), '--st-path', str(work / 'query/matrix.mtx.gz'),
            '--coordinates-path', str(work / 'coordinates.tsv'), '--st-cell-type-path', str(work / 'query_labels.tsv'),
            '--output-folder', str(work / 'upstream_output'), '--number-of-selected-spots', '10000',
            '--number-of-processors', '2', '--plot-off']
        env = {**os.environ, 'OPENBLAS_NUM_THREADS': '1', 'OMP_NUM_THREADS': '1', 'MKL_NUM_THREADS': '1',
            'NUMBA_CACHE_DIR': '/tmp/cytospace-numba', 'MPLCONFIGDIR': '/tmp/heart-mpl', 'PYTHONUNBUFFERED': '1'}
        start = time.time()
        print(f'RUN {stage}: {len(qi)} author-ventricular query segments, {len(ri)} independent reference nuclei; raw/default seed1/downsampling/10000 partitions', flush=True)
        with open(folder / 'upstream_stdout.txt', 'w') as log:
            subprocess.run(command, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
        assigned = pd.read_csv(work / 'upstream_output/assigned_locations.csv', float_precision='round_trip')
        assert assigned.SpotID.is_unique and set(assigned.SpotID) == set(labels.index[qi])
        assert assigned.CellType.eq('Ventricular_CM').all()
        qindex = labels.index.get_indexer(assigned.SpotID)
        rindex = metadata.index.get_indexer(assigned.OriginalCID)
        assert (rindex >= 0).all() and metadata.iloc[rindex].age.eq(stage).all()
        np.testing.assert_array_equal(assigned[['x_centroid', 'y_centroid']].to_numpy(), labels.iloc[qindex][['x_centroid', 'y_centroid']].to_numpy())
        assigned.to_csv(folder / 'assigned_locations.csv.gz', index=False, compression={'method': 'gzip', 'compresslevel': 3})
        out = assigned.copy()
        out['published_author_state'] = labels.iloc[qindex].author_atlas_cell_state.to_numpy()
        for key in ['sample', 'ventricular_group', 'mapping_family']:
            out['reference_' + key] = metadata.iloc[rindex][key].to_numpy()
        fit, shuffled_fit = [], []
        shuffled = np.random.default_rng(149).permutation(rindex)
        for offset in range(0, len(out), 4000):
            selected = slice(offset, offset + 4000)
            q = lognorm(query[qindex[selected]])
            r = lognorm(raw[rindex[selected]])
            s = lognorm(raw[shuffled[selected]])
            fit.extend(paired_r(q, r))
            shuffled_fit.extend(paired_r(q, s))
        out['original_reference_profile_pearson_r'] = fit
        out['within_ventricular_shuffled_profile_r'] = shuffled_fit
        out['reference_shared_molecules_gt1500'] = raw[rindex].sum(1) > 1500
        out.to_csv(folder / 'annotated_assignments.tsv.gz', sep='\t', index=False, na_rep='NA', compression={'method': 'gzip', 'compresslevel': 3})
        shutil.copyfile(work / 'upstream_output/log.txt', folder / 'upstream_log.txt')
        manifest = {'stage': stage, 'command': command, 'query_segments': len(qi), 'reference_nuclei': len(ri),
            'gene_count': len(genes), 'elapsed_seconds': time.time() - start,
            'query_cohort': 'Five ventricular-associated final states in public author WebAtlas, excluding only zero shared RNA',
            'source_of_reference': 'Our independent QC-reviewed snRNA reference, not original author reference objects',
            'source_of_spatial_labels': 'Retrieved published final author annotation; these are not verified original input labels',
            'subtype_labels_supplied': False, 'all_input_labels': 'Ventricular_CM',
            'default_parameters': {'seed': 1, 'downsample_off': False, 'max_shared_transcripts': 1500,
                'solver': 'lapjv', 'metric': 'Pearson_correlation', 'sampling': 'duplicates', 'partition_max': 10000},
            'input_order': 'Reference and query randomized with numpy default_rng1 before running upstream',
            'profile_fit_note': 'Diagnostic uses original raw reference, not the internally molecule-downsampled version; >1500 reference rows explicitly flagged',
            'source_query_sha256': digest(SOURCE / stage / f'{stage}_xenium_counts_and_candidate_mapping.h5ad'),
            'source_reference_sha256': digest(SOURCE / 'reference_shared_panel_counts.npz'),
            'assignments_sha256': digest(folder / 'assigned_locations.csv.gz'),
            'exact_author_pipeline_reproduction': False,
            'author_input_reference_barcodes_available': False,
            'parenthesized_comparison_warning': 'Changes include cohort, counts, seed and partitioning; does not isolate causal effects of any one setting'}
        (folder / 'run_manifest.json').write_text(json.dumps(manifest, indent=2))
        shutil.rmtree(work)
        print(f'DONE {stage}: {time.time()-start:.1f}s', flush=True)
        del data, query


if __name__ == '__main__':
    main()
