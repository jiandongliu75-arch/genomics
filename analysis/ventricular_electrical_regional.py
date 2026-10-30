"""Describe literature-guided electrical markers in existing ventricular groups.

No refitting, cell removal, imputation, or anatomical/functional assignment.
"""
from pathlib import Path
import argparse
import json

import h5py
import numpy as np
import pandas as pd
from anndata.io import read_elem, sparse_dataset
from scipy import sparse

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/ventricular_electrical_regional'
GENES = 'Slit2 Mpped2 Gja5 Etv1 Irx3 Cntn2 Scn10a Pcp4 Corin Mylk4 Tbx3 Irx5 Kcnd2 Kcnip2 Scn5a Gja1 Hey2'.split()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', default='/workspace/heart-atlas/results/cell_identities_subtypes/snrna_qc_annotated_subclustered.h5ad')
    args = parser.parse_args()
    OUT.mkdir(exist_ok=True)
    ann = pd.read_csv(ROOT / 'results/ventricular_microenvironment/cell_annotations_and_scores.csv.gz', index_col=0)
    with h5py.File(args.input) as handle:
        obs, var = read_elem(handle['obs']), read_elem(handle['var'])
        rows = obs.index.get_indexer(ann.index)
        assert (rows >= 0).all()
        names = var.gene_symbol.astype(str).to_numpy()
        print('Reading corrected counts', flush=True)
        matrix = sparse_dataset(handle['X'])[rows, :].tocsr().astype(np.float32)
        totals = np.asarray(matrix.sum(axis=1)).ravel()
        abundance = np.asarray(matrix.sum(axis=0)).ravel()
        selected, mapping = [], []
        for gene in GENES:
            matches = np.flatnonzero(names == gene)
            if not len(matches):
                raise ValueError(f'Missing gene: {gene}')
            feature = int(matches[np.argmax(abundance[matches])])
            selected.append(feature)
            mapping.append({'gene': gene, 'gene_id': var.index[feature], 'matching_features': ';'.join(var.index[matches])})
        corrected = matrix[:, selected].toarray()
        del matrix
        print('Reading original counts', flush=True)
        matrix = sparse_dataset(handle['layers/counts_raw'])[rows, :].tocsr()
        rawtotals = np.asarray(matrix.sum(axis=1)).ravel()
        raw = matrix[:, selected].toarray()
        del matrix
    log = np.log1p(corrected * 10000 / totals[:, None])
    rawlog = np.log1p(raw * 10000 / rawtotals[:, None])
    residual = log.copy()
    rawresidual = rawlog.copy()
    for idx in ann.groupby(['age', 'sample'], observed=True).indices.values():
        sub = ann.iloc[idx]
        cov = np.column_stack([np.log1p(sub.raw_total_counts), np.log1p(sub.raw_n_genes_by_counts), sub.decontx_contamination, sub.raw_pct_counts_mt, sub.cycling_marker_count])
        cov = (cov - cov.mean(axis=0)) / np.maximum(cov.std(axis=0), 1e-8)
        residual[idx] -= cov @ np.linalg.lstsq(cov, log[idx] - log[idx].mean(axis=0), rcond=None)[0]
        rawresidual[idx] -= cov @ np.linalg.lstsq(cov, rawlog[idx] - rawlog[idx].mean(axis=0), rcond=None)[0]
    profiles = []
    for fields in [['age'], ['age', 'sample'], ['age', 'microenvironment_cluster'], ['age', 'sample', 'microenvironment_cluster']]:
        for key, idx in ann.groupby(fields, observed=True).indices.items():
            key = key if isinstance(key, tuple) else (key,)
            meta = dict(zip(fields, key))
            for j, gene in enumerate(GENES):
                profiles.append({**meta, 'grouping': '+'.join(fields), 'gene': gene, 'nuclei': len(idx), 'mean_log1p_10k': float(log[idx, j].mean()), 'mean_raw_log1p_10k': float(rawlog[idx, j].mean()), 'fraction_corrected_ge1': float((corrected[idx, j] >= 1).mean()), 'fraction_raw_ge1': float((raw[idx, j] >= 1).mean())})
    pd.DataFrame(profiles).to_csv(OUT / 'gene_profiles.tsv', sep='\t', index=False, na_rep='NA')
    contrasts, library_rows = [], []
    for population in ['expanded', 'original_core_only']:
        selected_population = np.ones(len(ann), dtype=bool) if population == 'expanded' else ann.gate_origin.eq('core_ventricular').to_numpy()
        for (age, cluster), allidx in ann.groupby(['age', 'microenvironment_cluster'], observed=True).indices.items():
            clusteridx = allidx[selected_population[allidx]]
            effects = []
            for sample in sorted(ann.loc[ann.age.eq(age), 'sample'].unique()):
                libmask = ann.age.eq(age).to_numpy() & ann['sample'].eq(sample).to_numpy() & selected_population
                ingroup = np.flatnonzero(libmask & ann.microenvironment_cluster.eq(cluster).to_numpy())
                rest = np.flatnonzero(libmask & ~ann.microenvironment_cluster.eq(cluster).to_numpy())
                eligible = len(ingroup) >= 30 and len(rest) >= 30
                delta = log[ingroup].mean(axis=0) - log[rest].mean(axis=0) if eligible else np.full(len(GENES), np.nan)
                rawdelta = rawlog[ingroup].mean(axis=0) - rawlog[rest].mean(axis=0) if eligible else np.full(len(GENES), np.nan)
                adjusted = residual[ingroup].mean(axis=0) - residual[rest].mean(axis=0) if eligible else np.full(len(GENES), np.nan)
                rawadjusted = rawresidual[ingroup].mean(axis=0) - rawresidual[rest].mean(axis=0) if eligible else np.full(len(GENES), np.nan)
                if eligible:
                    effects.append((delta, rawdelta, adjusted, rawadjusted))
                for j, gene in enumerate(GENES):
                    library_rows.append({'population': population, 'age': age, 'sample': sample, 'cluster': cluster, 'gene': gene, 'nuclei': len(ingroup), 'other_nuclei': len(rest), 'eligible': eligible, 'corrected_delta': delta[j], 'raw_delta': rawdelta[j], 'QC_cycle_adjusted_delta': adjusted[j], 'raw_QC_cycle_adjusted_delta': rawadjusted[j]})
            effects = np.asarray(effects)
            for j, gene in enumerate(GENES):
                two = len(effects) == 2
                direction = np.sign(effects[:, 0, j]) if len(effects) else np.array([])
                contrasts.append({'population': population, 'age': age, 'cluster': cluster, 'gene': gene, 'nuclei': len(clusteridx), 'reference_nuclei': int(ann.iloc[clusteridx].gate_origin.eq('Slit2_overlap_reference').sum()), 'eligible_libraries': len(effects), 'corrected_delta': float(effects[:, 0, j].mean()) if len(effects) else np.nan, 'raw_delta': float(effects[:, 1, j].mean()) if len(effects) else np.nan, 'QC_cycle_adjusted_delta': float(effects[:, 2, j].mean()) if len(effects) else np.nan, 'raw_QC_cycle_adjusted_delta': float(effects[:, 3, j].mean()) if len(effects) else np.nan, 'both_libraries_corrected_positive': bool(two and (direction > 0).all()), 'both_libraries_corrected_negative': bool(two and (direction < 0).all()), 'both_libraries_raw_same_direction': bool(two and (direction == np.sign(effects[:, 1, j])).all() and abs(direction.sum()) == 2), 'both_libraries_QC_adjusted_same_direction': bool(two and (direction == np.sign(effects[:, 2, j])).all() and abs(direction.sum()) == 2)})
    combined = pd.DataFrame(contrasts)
    combined[combined.population.eq('expanded')].to_csv(OUT / 'cluster_gene_contrasts.tsv', sep='\t', index=False, na_rep='NA')
    combined[combined.population.eq('original_core_only')].to_csv(OUT / 'core_only_cluster_gene_contrasts.tsv', sep='\t', index=False, na_rep='NA')
    pd.DataFrame(library_rows).to_csv(OUT / 'library_cluster_gene_contrasts.tsv', sep='\t', index=False, na_rep='NA')
    pd.DataFrame(mapping).to_csv(OUT / 'gene_mapping.tsv', sep='\t', index=False, na_rep='NA')
    columns = ['age', 'sample', 'microenvironment_cluster', 'gate_origin']
    ann[columns].to_csv(OUT / 'cell_metadata.tsv.gz', sep='\t')
    expression = pd.DataFrame(log, index=ann.index, columns=GENES)
    expression.to_csv(OUT / 'targeted_log_expression.csv.gz', float_format='%.7g')
    sparse.save_npz(OUT / 'targeted_corrected_counts.npz', sparse.csr_matrix(corrected))
    sparse.save_npz(OUT / 'targeted_original_counts.npz', sparse.csr_matrix(raw))
    pd.DataFrame({'corrected_whole_nucleus_total': totals, 'original_whole_nucleus_total': rawtotals}, index=ann.index).to_csv(OUT / 'normalization_totals.tsv.gz', sep='\t')
    (OUT / 'methods.json').write_text(json.dumps({'input': args.input, 'nuclei': len(ann), 'genes': GENES, 'clusters': 28, 'expression': 'log1p(count / whole-nucleus count total * 10000), independently for fractional corrected and original integer counts.', 'contrast': 'Cluster minus all other retained CMs in the same age and original library; require >=30 nuclei in each group. Summary is equal-weight mean over eligible libraries; both-library flags require two eligible libraries. No p-values or biological-replicate inference.', 'QC_sensitivity': 'Within-library expression residuals after linear removal of log total RNA, log detected genes, contamination, mitochondrial fraction and cycling-marker count. May remove genuine biological mediators; conditional sensitivity only.', 'core_sensitivity': 'Repeat contrasts using only original core ventricular nuclei, preserving original clusters. Some reference-rich clusters become ineligible.', 'selection_limits': 'Existing groups were interpreted using Slit2/Mpped2/Gja5/Etv1 panel and historical conduction-associated gates; these genes and related markers are not independent validation. No cell identity changed.', 'matrix_files': 'Sparse NPZ count matrices: rows in cell_metadata.tsv.gz order, columns in gene_mapping.tsv order. Only 17 selected features, not full transcriptome. Corrected counts remain fractional.'}, indent=2))
    print(f'DONE: {len(ann)} nuclei; {len(GENES)} genes; unchanged 28 groups.', flush=True)


if __name__ == '__main__':
    main()
