"""Run unmodified upstream CytoSPACE in single-cell mode on matched heart stages.

Use the main heart-atlas Python for preparation; subprocesses use a separate,
version-pinned Python 3.9 environment for the upstream implementation.
"""
from pathlib import Path
import argparse
import gzip
import hashlib
import importlib.metadata
import json
import os
import shutil
import subprocess
import time

import anndata as ad
import numpy as np
import pandas as pd
from scipy import io, sparse

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'results/xenium_ventricular_mapping'
OUT = ROOT / 'results/cytospace_xenium'
BASE = Path('/workspace/heart-atlas')
CACHE = BASE / 'results/cytospace_cache'
CYTO = BASE / 'venv-cytospace/bin/cytospace'
STAGES = ['P0', 'P7', 'P14', 'P21']
RUNS = ['corrected_seed42', 'corrected_seed43', 'corrected_seed44',
        'original_seed42', 'library1_seed42', 'library2_seed42', 'heldout_seed42']
LINEAGES = {'Ventricular-enriched CM': 'Ventricular_CM',
            'Atrial-marker-enriched CM': 'Atrial_CM', 'Endothelial': 'Endothelial',
            'Fibroblast': 'Fibroblast', 'Pericyte': 'Pericyte',
            'Smooth muscle': 'Smooth_muscle', 'Epicardial': 'Epicardial',
            'Immune': 'Immune'}
REF_LINEAGES = {**LINEAGES, 'Atrial-like CM': 'Atrial_CM'}
# These genes are not positive/exclusion markers in the ventricular gate.
# Total-count/QC selection is fixed. This tests recovery conditional on the
# supplied broad identities, not end-to-end blinded annotation accuracy.
HOLDOUT = 'Slit2 Ankrd1 Cntn2 Gja5 Ddc Mki67 Top2a Sema3a Ephb1 Thsd4 Cpne5 Pola1 Rad51b Cd36 Cpt1b Ppara Ppargc1a Ppargc1b Kcnip2 Scn5a Hey2'.split()


def digest(path):
    h = hashlib.sha256()
    with open(path, 'rb') as handle:
        for part in iter(lambda: handle.read(8 << 20), b''):
            h.update(part)
    return h.hexdigest()


def write_matrix(folder, matrix, genes, ids):
    folder.mkdir(parents=True, exist_ok=True)
    with gzip.open(folder / 'matrix.mtx.gz', 'wb', compresslevel=2) as handle:
        io.mmwrite(handle, sparse.coo_matrix(matrix.T), symmetry='general')
    pd.Series(genes).to_csv(folder / 'genes.tsv', sep='\t', index=False, header=False)
    pd.Series(ids).to_csv(folder / 'cells.tsv', sep='\t', index=False, header=False)


def lognorm(matrix):
    matrix = np.asarray(matrix, dtype=np.float64)
    return np.log2(1 + 1e6 * matrix / matrix.sum(axis=1, keepdims=True))


def paired_r(first, second):
    first = first - first.mean(axis=1, keepdims=True)
    second = second - second.mean(axis=1, keepdims=True)
    denom = np.sqrt((first * first).sum(1) * (second * second).sum(1))
    return np.divide((first * second).sum(1), denom, out=np.zeros(len(first)), where=denom > 0)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--stages', nargs='+', default=STAGES, choices=STAGES)
    parser.add_argument('--runs', nargs='+', default=RUNS, choices=RUNS)
    parser.add_argument('--chunk-size', type=int, default=3000)
    parser.add_argument('--processors', type=int, default=4)
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    CACHE.mkdir(parents=True, exist_ok=True)
    metadata = pd.read_csv(SOURCE / 'reference_metadata.tsv.gz', sep='\t', index_col=0)
    bundle = np.load(SOURCE / 'reference_shared_panel_counts.npz')
    genes = bundle['genes'].astype(str)
    assert np.array_equal(bundle['cell_ids'], metadata.index.to_numpy())
    corrected, original = bundle['corrected'], bundle['original']
    heldout = [gene for gene in HOLDOUT if gene in genes]
    source_sha = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=BASE / 'tools/cytospace', text=True).strip()
    assert source_sha == 'bb7fe12ab4466ef5e38dcac4c3ece07012060d60', 'Use the pinned upstream commit'
    subprocess.run(['git', 'diff', '--exit-code', 'HEAD', '--', 'cytospace'], cwd=BASE / 'tools/cytospace', check=True, stdout=subprocess.DEVNULL)
    for stage in args.stages:
        stage_out = OUT / stage
        stage_out.mkdir(exist_ok=True)
        print(f'LOAD {stage}', flush=True)
        data = ad.read_h5ad(SOURCE / stage / f'{stage}_xenium_counts_and_candidate_mapping.h5ad')
        obs = data.obs.copy()
        broad = obs.spatial_lineage_candidate.astype(str).map(LINEAGES)
        eligible = obs.pass_expression_QC & obs.pass_transfer_information & broad.notna()
        obs['cytospace_input_lineage'] = broad.fillna('Not included')
        obs['cytospace_eligible'] = eligible
        obs['cytospace_exclusion_reason'] = np.where(eligible, 'Included', np.where(
            ~(obs.pass_expression_QC & obs.pass_transfer_information), 'QC/information threshold',
            'Unresolved chamber/lineage or CM/EC overlap review'))
        obs.to_csv(stage_out / 'all_cells_eligibility.tsv.gz', sep='\t', na_rep='NA')
        query = data[:, genes].X.toarray().astype(np.float32)
        ref_eligible = metadata.age.eq(stage) & metadata.broad_class.isin(REF_LINEAGES)
        ref_eligible &= (corrected.sum(1) >= 20) & ((corrected > 0).sum(1) >= 8)
        ref_eligible &= (original.sum(1) >= 20) & ((original > 0).sum(1) >= 8)
        libraries = sorted(metadata.loc[ref_eligible, 'sample'].unique())
        assert len(libraries) == 2
        for run in args.runs:
            destination = stage_out / run
            if (destination / 'run_manifest.json').exists():
                print(f'SKIP_COMPLETED {stage} {run}', flush=True)
                continue
            destination.mkdir(exist_ok=True)
            seed = int(run.rsplit('seed', 1)[1])
            selected_ref = np.flatnonzero(ref_eligible)
            if run.startswith('library'):
                library = libraries[int(run[7]) - 1]
                selected_ref = selected_ref[metadata.iloc[selected_ref]['sample'].eq(library)]
            selected_query = np.flatnonzero(eligible)
            # The upstream labelled single-cell mode does not shuffle partition
            # order itself. Randomize both matrices here, preserving ID alignment.
            rng = np.random.default_rng(seed)
            selected_ref = rng.permutation(selected_ref)
            selected_query = rng.permutation(selected_query)
            ref = metadata.iloc[selected_ref].copy()
            ref['cytospace_input_lineage'] = ref.broad_class.map(REF_LINEAGES)
            included_types = set(broad.iloc[selected_query]) & set(ref.cytospace_input_lineage)
            selected_query = selected_query[broad.iloc[selected_query].isin(included_types)]
            selected_ref = selected_ref[ref.cytospace_input_lineage.isin(included_types)]
            ref = metadata.iloc[selected_ref].copy()
            ref['cytospace_input_lineage'] = ref.broad_class.map(REF_LINEAGES)
            use = np.array([not (run.startswith('heldout') and gene in heldout) for gene in genes])
            train_genes = genes[use]
            ref_counts = original if run.startswith('original') else corrected
            work = CACHE / stage / run
            work.mkdir(parents=True, exist_ok=True)
            write_matrix(work / 'reference', ref_counts[selected_ref][:, use], train_genes, ref.index)
            write_matrix(work / 'query', query[selected_query][:, use], train_genes, obs.index[selected_query])
            ref[['cytospace_input_lineage']].to_csv(work / 'reference_labels.tsv', sep='\t', index_label='CellID')
            pd.DataFrame({'CellType': broad.iloc[selected_query].to_numpy()}, index=obs.index[selected_query]).to_csv(work / 'query_labels.tsv', sep='\t', index_label='SpotID')
            obs.iloc[selected_query][['x_centroid', 'y_centroid']].to_csv(work / 'coordinates.tsv', sep='\t', index_label='SpotID')
            ref.to_csv(destination / 'reference_pool.tsv.gz', sep='\t', na_rep='NA')
            pd.Series(obs.index[selected_query], name='SpotID').to_csv(destination / 'query_input_order.tsv.gz', sep='\t', index=False)
            command = [str(CYTO), '--single-cell', '--scRNA-path', str(work / 'reference/matrix.mtx.gz'),
                       '--cell-type-path', str(work / 'reference_labels.tsv'), '--st-path', str(work / 'query/matrix.mtx.gz'),
                       '--coordinates-path', str(work / 'coordinates.tsv'), '--st-cell-type-path', str(work / 'query_labels.tsv'),
                       '--output-folder', str(work / 'upstream_output'), '--number-of-selected-spots', str(args.chunk_size),
                       '--number-of-processors', str(args.processors), '--solver-method', 'lapjv',
                       '--distance-metric', 'Pearson_correlation', '--sampling-method', 'duplicates',
                       '--seed', str(seed), '--downsample-off', '--plot-off']
            env = {**os.environ, 'OPENBLAS_NUM_THREADS': '1', 'OMP_NUM_THREADS': '1', 'MKL_NUM_THREADS': '1',
                   'NUMBA_CACHE_DIR': '/tmp/cytospace-numba', 'MPLCONFIGDIR': '/tmp/heart-mpl', 'PYTHONUNBUFFERED': '1'}
            start = time.time()
            print(f'RUN {stage} {run}: {len(selected_query)} query / {len(ref)} reference / {len(train_genes)} genes', flush=True)
            with open(destination / 'upstream_stdout.txt', 'w') as log:
                subprocess.run(command, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
            assigned = pd.read_csv(work / 'upstream_output/assigned_locations.csv', float_precision='round_trip')
            assert assigned.SpotID.is_unique and len(assigned) == len(selected_query)
            assert set(assigned.SpotID) == set(obs.index[selected_query])
            qi = obs.index.get_indexer(assigned.SpotID)
            ri = metadata.index.get_indexer(assigned.OriginalCID)
            assert (ri >= 0).all()
            assert np.array_equal(assigned[['x_centroid', 'y_centroid']].to_numpy(), obs.iloc[qi][['x_centroid', 'y_centroid']].to_numpy())
            assert np.array_equal(assigned.CellType.to_numpy(), broad.iloc[qi].to_numpy())
            assigned.to_csv(destination / 'assigned_locations.csv.gz', index=False)
            annotations = assigned.copy()
            for field in ['sample', 'broad_class', 'ventricular_group', 'mapping_family']:
                annotations['reference_' + field] = metadata.iloc[ri][field].to_numpy()
            # These are expression-fit diagnostics, not confidence probabilities.
            correlations, shuffled = [], []
            shuffle_ri = ri.copy()
            shuffle_rng = np.random.default_rng(8491)
            for lineage in assigned.CellType.unique():
                indices = np.flatnonzero(assigned.CellType.eq(lineage))
                shuffle_ri[indices] = shuffle_rng.permutation(ri[indices])
            for offset in range(0, len(assigned), 4000):
                sl = slice(offset, offset + 4000)
                qq = lognorm(query[qi[sl]][:, use])
                rr = lognorm(ref_counts[ri[sl]][:, use])
                ss = lognorm(ref_counts[shuffle_ri[sl]][:, use])
                correlations.extend(paired_r(qq, rr))
                shuffled.extend(paired_r(qq, ss))
            annotations['mapping_gene_pearson_r'] = correlations
            annotations['within_lineage_shuffled_pearson_r'] = shuffled
            annotations.to_csv(destination / 'annotated_assignments.tsv.gz', sep='\t', index=False, na_rep='NA')
            shutil.copyfile(work / 'upstream_output/log.txt', destination / 'upstream_log.txt')
            manifest = {'stage': stage, 'run': run, 'seed': seed, 'upstream_version': '1.1.0', 'upstream_commit': source_sha,
                        'upstream_source_modified': False, 'command': command, 'elapsed_seconds_including_checks': time.time() - start,
                        'mapping_genes': train_genes.tolist(), 'withheld_genes': heldout if run.startswith('heldout') else [],
                        'query_cells': len(assigned), 'reference_pool_nuclei': len(ref), 'reference_libraries': sorted(ref['sample'].unique()),
                        'reference_counts': 'original integer counts' if run.startswith('original') else 'DecontX-corrected unnormalized fractional count estimates',
                        'reference_sampling': 'upstream duplicates: every available nucleus once, then sample with replacement when more query cells are needed; no subtype labels/fractions supplied',
                        'broad_labels': 'fixed input lineage/chamber candidates from preceding mapping; not independently re-identified by CytoSPACE',
                        'partitioning': 'per broad lineage, randomized query/reference input order, then upstream chunks; coordinates never enter expression costs',
                        'chunk_size': args.chunk_size, 'downsampling': 'off; normalization is upstream log2(1+1e6*counts/shared-panel total)',
                        'assigned_expression': 'upstream redundant shared-panel copies omitted from publication; recover each column by OriginalCID using the published reference panel',
                        'source_query_sha256': digest(SOURCE / stage / f'{stage}_xenium_counts_and_candidate_mapping.h5ad'),
                        'source_reference_sha256': digest(SOURCE / 'reference_shared_panel_counts.npz'),
                        'assignments_sha256': digest(destination / 'assigned_locations.csv.gz')}
            (destination / 'run_manifest.json').write_text(json.dumps(manifest, indent=2))
            # These are our transient reproducible inputs and redundant upstream
            # reference expression copies, not original data. Bound disk usage.
            shutil.rmtree(work)
            print(f'DONE {stage} {run}: {time.time()-start:.1f}s', flush=True)
        del data, query


if __name__ == '__main__':
    main()
