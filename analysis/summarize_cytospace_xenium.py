"""Audit forced CytoSPACE matches, filter stable candidates, and test held-out genes."""
from pathlib import Path
import argparse
import json
import time

import anndata as ad
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, pearsonr
from sklearn.metrics import roc_auc_score

from run_cytospace_xenium import SOURCE, OUT, STAGES, RUNS, lognorm


def r_safe(x, y, rank=False):
    if len(x) < 3 or np.std(x) == 0 or np.std(y) == 0:
        return np.nan
    return float((spearmanr if rank else pearsonr)(x, y).statistic)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--stages', nargs='+', default=STAGES, choices=STAGES)
    parser.add_argument('--wait-for-runs', action='store_true', help='Wait for ongoing stage runs, up to 30 minutes per stage')
    args = parser.parse_args()
    ref = pd.read_csv(SOURCE / 'reference_metadata.tsv.gz', sep='\t', index_col=0)
    panel = np.load(SOURCE / 'reference_shared_panel_counts.npz')
    genes = panel['genes'].astype(str)
    reference_counts = panel['corrected']
    catalog = pd.read_csv(SOURCE / 'reference_group_catalog.tsv', sep='\t').set_index('microenvironment_cluster')
    recoverable = pd.read_csv(OUT / 'reference_panel_benchmark/panel_recoverability.tsv', sep='\t')
    marker_values = pd.read_csv(OUT / 'snRNA_reference_marker_values_NOT_Xenium_measurements.tsv.gz', sep='\t', index_col=0)
    summaries, counts, comparisons, fits, priors, heldout_metrics, grid_values = [], [], [], [], [], [], []
    group_markers = []
    for stage in args.stages:
        folder = OUT / stage
        if args.wait_for_runs:
            start = time.monotonic()
            missing = [run for run in RUNS if not (folder / run / 'run_manifest.json').exists()]
            if missing:
                print(f'WAIT {stage}: {missing}', flush=True)
            while missing:
                if time.monotonic() - start > 1800:
                    raise TimeoutError(f'Unfinished {stage} runs: {missing}')
                time.sleep(5)
                missing = [run for run in RUNS if not (folder / run / 'run_manifest.json').exists()]
        print(f'SUMMARIZE {stage}', flush=True)
        obs = pd.read_csv(folder / 'all_cells_eligibility.tsv.gz', sep='\t', index_col=0, float_precision='round_trip')
        data = ad.read_h5ad(SOURCE / stage / f'{stage}_xenium_counts_and_candidate_mapping.h5ad')
        assert np.array_equal(data.obs_names, obs.index)
        query_counts = data[:, genes].X.toarray().astype(np.float32)
        runs = {run: pd.read_csv(folder / run / 'annotated_assignments.tsv.gz', sep='\t', float_precision='round_trip').set_index('SpotID') for run in RUNS}
        primary = runs['corrected_seed42']
        selected = obs.index[obs.cytospace_eligible]
        assert all(set(df.index) == set(selected) for df in runs.values())
        vent = obs.spatial_lineage_candidate.eq('Ventricular-enriched CM')
        keep = ['x_centroid', 'y_centroid', 'stage', 'section', 'transcript_counts', 'cell_area', 'nucleus_area',
                'shared_panel_transcripts', 'shared_panel_detected_genes', 'missing_nucleus_area_flag',
                'pass_expression_QC', 'pass_transfer_information', 'CM_EC_program_overlap_review_flag',
                'spatial_lineage_candidate', 'spatial_mapped_family', 'spatial_mapped_fine_group',
                'measured_Ankrd1_program_positive', 'measured_Slit2_program_positive',
                'measured_conduction_overlap_positive', 'measured_cycling_program_positive',
                'distance_to_Npr3_endothelial_landmark_um', 'distance_to_Upk3b_Wt1_surface_landmark_um',
                'cytospace_eligible', 'cytospace_input_lineage', 'cytospace_exclusion_reason']
        out = obs[keep].copy()
        for run, df in runs.items():
            for source, suffix in [('OriginalCID', 'reference_cell_id'), ('reference_sample', 'reference_library'),
                                   ('reference_mapping_family', 'family'), ('reference_ventricular_group', 'group'),
                                   ('mapping_gene_pearson_r', 'pearson_r')]:
                out[f'cyto_{run}_{suffix}'] = df[source].reindex(obs.index)
            for lineage, indices in df.groupby('CellType').groups.items():
                fits.append({'stage': stage, 'run': run, 'input_lineage': lineage, 'cells': len(indices),
                    'median_assigned_pearson_r': df.loc[indices, 'mapping_gene_pearson_r'].median(),
                    'median_within_lineage_shuffled_r': df.loc[indices, 'within_lineage_shuffled_pearson_r'].median()})
        # Stability is not a calibrated confidence score.
        for resolution in ['family', 'group']:
            columns = [f'cyto_corrected_seed{seed}_{resolution}' for seed in [42, 43, 44]]
            values = out.loc[vent, columns].fillna('NA')
            mode = values.mode(axis=1)[0]
            n_agree = values.eq(mode, axis=0).sum(axis=1)
            out.loc[vent, f'cyto_seed_consensus_{resolution}'] = mode
            out.loc[vent, f'cyto_seed_consensus_{resolution}_votes'] = n_agree
            out.loc[vent, f'cyto_primary_matches_consensus_{resolution}'] = values.iloc[:, 0].eq(mode)
            for name in ['original_seed42', 'library1_seed42', 'library2_seed42']:
                out.loc[vent, f'cyto_{name}_{resolution}_agreement'] = out.loc[vent, f'cyto_{name}_{resolution}'].eq(mode)
        stage_benchmark = recoverable[recoverable.stage.eq(stage)]
        good_families = set(stage_benchmark.loc[stage_benchmark.resolution.eq('mapping_family') & stage_benchmark.panel_recoverable, 'label'])
        good_groups = set(stage_benchmark.loc[stage_benchmark.resolution.eq('ventricular_group') & stage_benchmark.panel_recoverable, 'label'])
        family = out.cyto_seed_consensus_family
        group = out.cyto_seed_consensus_group
        family_stable = vent & (out.cyto_seed_consensus_family_votes >= 2) & out.cyto_primary_matches_consensus_family.eq(True)
        group_stable = vent & (out.cyto_seed_consensus_group_votes >= 2) & out.cyto_primary_matches_consensus_group.eq(True)
        for run in ['original_seed42', 'library1_seed42', 'library2_seed42']:
            family_stable &= out[f'cyto_{run}_family_agreement'].eq(True)
            group_stable &= out[f'cyto_{run}_group_agreement'].eq(True)
        family_stable &= (out.cyto_corrected_seed42_pearson_r >= .30) & (obs.shared_panel_transcripts >= 50)
        marker_support = (~family.eq('Cycling-associated') | (obs.measured_cycling_genes >= 2))
        marker_support &= (~family.eq('Conduction-program candidate') | (obs.measured_conduction_overlap_genes >= 2))
        family_ok = family_stable & family.isin(good_families) & marker_support
        query_marker_n = pd.Series(0, index=obs.index)
        for label in catalog.index[catalog.stage.eq(stage)]:
            shared = [g for g in str(catalog.loc[label, 'both_library_markers']).split(';') if g in genes]
            ci = np.flatnonzero(vent & group.eq(label))
            indices = [np.flatnonzero(genes == g)[0] for g in shared]
            if len(ci) and indices:
                query_marker_n.iloc[ci] = (query_counts[ci][:, indices] > 0).sum(1)
            group_markers.append({'stage': stage, 'group': label, 'measured_both_library_marker_genes': ';'.join(shared), 'n_measured': len(shared)})
        group_ok = group_stable & group.isin(good_groups) & family_ok & (obs.shared_panel_transcripts >= 100)
        group_ok &= (query_marker_n >= 3) & group.map(catalog.mapping_family).eq(family)
        out['cyto_robustness_only_family_pass'] = family_stable & marker_support
        out['cyto_robustness_only_group_pass'] = (group_stable & marker_support &
            (obs.shared_panel_transcripts >= 100) & (query_marker_n >= 3))
        out['cyto_panel_recoverable_family'] = family.isin(good_families)
        out['cyto_panel_recoverable_group'] = group.isin(good_groups)
        out['cyto_consensus_group_marker_genes_detected'] = query_marker_n
        out['cyto_stable_family_candidate'] = np.where(family_ok, family, np.where(vent, 'Not supported by all checks', 'Not a ventricular input'))
        out['cyto_stable_group_candidate'] = np.where(group_ok, group, np.where(vent, 'Not supported by all checks', 'Not a ventricular input'))
        # Gene values carried by the assigned snRNA nuclei. All column names
        # explicitly say prediction; these are never merged into measured counts.
        for gene in ['Kcnd2', 'Irx5', 'Irx3', 'Etv1', 'Mpped2', 'Kcnip2']:
            predictions = []
            for seed in [42, 43, 44]:
                ids = runs[f'corrected_seed{seed}'].loc[out.index[vent], 'OriginalCID']
                predictions.append(marker_values.loc[ids, f'{gene}_corrected_log1pCP10k'].to_numpy())
            predictions = np.array(predictions).T
            out.loc[vent, f'predicted_reference_{gene}_mean_log1pCP10k'] = predictions.mean(1)
            out.loc[vent, f'predicted_reference_{gene}_range_across_3seeds'] = np.ptp(predictions, axis=1)
        out.to_csv(folder / 'cytospace_cell_annotations.tsv.gz', sep='\t', na_rep='NA')
        for field in ['cyto_corrected_seed42_group', 'cyto_seed_consensus_group', 'cyto_stable_group_candidate',
                      'cyto_corrected_seed42_family', 'cyto_stable_family_candidate']:
            for label, n in out.loc[vent, field].value_counts(dropna=False).items():
                counts.append({'stage': stage, 'assignment_field': field, 'label': label, 'cells': n})
        for run, df in runs.items():
            labels = df.loc[df.CellType.eq('Ventricular_CM'), 'reference_ventricular_group'].value_counts()
            pool = pd.read_csv(folder / run / 'reference_pool.tsv.gz', sep='\t', index_col=0)
            pool = pool.loc[pool.broad_class.eq('Ventricular-enriched CM'), 'ventricular_group'].value_counts()
            for label in sorted(set(labels.index) | set(pool.index)):
                priors.append({'stage': stage, 'run': run, 'group': label, 'reference_pool_nuclei': pool.get(label, 0),
                    'reference_pool_fraction': pool.get(label, 0) / pool.sum(), 'assigned_query_cells': labels.get(label, 0),
                    'assigned_query_fraction': labels.get(label, 0) / labels.sum()})
        old_family = obs.spatial_mapped_family
        old_fine = obs.spatial_mapped_fine_group
        previously_assigned = vent & ~old_family.eq('Unresolved ventricular state')
        both = previously_assigned & family_ok
        comparisons.append({'stage': stage, 'prior_family_assigned': int(previously_assigned.sum()),
            'prior_assigned_vs_native_cytospace_family_agreement': float(out.loc[previously_assigned, 'cyto_corrected_seed42_family'].eq(old_family[previously_assigned]).mean()),
            'both_methods_filtered_family_assigned': int(both.sum()),
            'filtered_family_agreement': float(out.loc[both, 'cyto_stable_family_candidate'].eq(old_family[both]).mean()) if both.any() else np.nan,
            'both_methods_filtered_group_assigned': int((group_ok & ~old_fine.eq('Not resolved')).sum())})
        pd.crosstab(old_family[vent], out.loc[vent, 'cyto_corrected_seed42_family']).to_csv(folder / 'prior_vs_cytospace_native_family.tsv', sep='\t')
        summaries.append({'stage': stage, 'original_segments': len(obs), 'cytospace_included_segments': int(obs.cytospace_eligible.sum()),
            'excluded_retained_segments': int((~obs.cytospace_eligible).sum()), 'ventricular_inputs': int(vent.sum()),
            'robustness_only_families': int((family_stable & marker_support).sum()),
            'panel_supported_stable_family_candidates': int(family_ok.sum()), 'panel_supported_stable_fine_candidates': int(group_ok.sum()),
            'ventricular_seed_family_majority_fraction': float((out.loc[vent, 'cyto_seed_consensus_family_votes'] >= 2).mean()),
            'ventricular_seed_group_majority_fraction': float((out.loc[vent, 'cyto_seed_consensus_group_votes'] >= 2).mean()),
            'ventricular_all3_seed_same_family_fraction': float(out.loc[vent, 'cyto_seed_consensus_family_votes'].eq(3).mean()),
            'ventricular_all3_seed_same_group_fraction': float(out.loc[vent, 'cyto_seed_consensus_group_votes'].eq(3).mean()),
            'ventricular_median_fit_r': float(primary.loc[primary.CellType.eq('Ventricular_CM'), 'mapping_gene_pearson_r'].median())})
        # Held-out gene recovery, exclusively within the fixed ventricular cohort.
        held = runs['heldout_seed42'].loc[out.index[vent]]
        qi = obs.index.get_indexer(held.index)
        ri = ref.index.get_indexer(held.OriginalCID)
        qnorm, rnorm = lognorm(query_counts[qi]), lognorm(reference_counts[ri])
        manifest = json.loads((folder / 'heldout_seed42/run_manifest.json').read_text())
        covariates = np.column_stack([np.ones(len(qi)), np.log1p(obs.iloc[qi].shared_panel_transcripts),
            np.log1p(obs.iloc[qi].cell_area), np.log1p(obs.iloc[qi].shared_panel_detected_genes),
            obs.iloc[qi].missing_nucleus_area_flag.astype(float)])
        cov_q, _ = np.linalg.qr(covariates, mode='reduced')
        coords = obs.iloc[qi][['x_centroid', 'y_centroid']].to_numpy()
        grid = np.floor(coords / 250).astype(int)
        bins = pd.DataFrame({'grid_x': grid[:, 0], 'grid_y': grid[:, 1]})
        for gene in manifest['withheld_genes']:
            j = np.flatnonzero(genes == gene)[0]
            x, y = rnorm[:, j], qnorm[:, j]
            positive = query_counts[qi, j] > 0
            x_res = x - cov_q @ (cov_q.T @ x)
            y_res = y - cov_q @ (cov_q.T @ y)
            binned = bins.assign(measured=y, predicted=x).groupby(['grid_x', 'grid_y']).agg(
                measured_mean=('measured', 'mean'), predicted_mean=('predicted', 'mean'), cells=('measured', 'size')).reset_index()
            binned = binned[binned.cells >= 20]
            binned['stage'], binned['gene'] = stage, gene
            grid_values.append(binned)
            shuffle_x = np.random.default_rng(491).permutation(x)
            heldout_metrics.append({'stage': stage, 'gene': gene, 'ventricular_cells': len(qi), 'measured_positive_cells': int(positive.sum()),
                'cell_spearman_rho': r_safe(x, y, rank=True), 'cell_log_normalized_pearson_r': r_safe(x, y),
                'depth_area_gene_count_adjusted_r': r_safe(x_res, y_res),
                'detection_AUC': roc_auc_score(positive, x) if len(np.unique(positive)) == 2 else np.nan,
                'within_ventricular_shuffle_rho': r_safe(shuffle_x, y, rank=True),
                'grid250um_bins_ge20cells': len(binned),
                'grid250um_spearman_rho': r_safe(binned.predicted_mean, binned.measured_mean, rank=True),
                'normalization': 'upstream log2(CP1M+1), calculated over full shared panel for evaluation only'})
        del data, query_counts
    pd.DataFrame(summaries).to_csv(OUT / 'stage_summary.tsv', sep='\t', index=False)
    pd.DataFrame(counts).to_csv(OUT / 'ventricular_assignment_counts.tsv', sep='\t', index=False)
    pd.DataFrame(comparisons).to_csv(OUT / 'comparison_with_previous_mapping.tsv', sep='\t', index=False)
    pd.DataFrame(fits).to_csv(OUT / 'expression_fit_diagnostics.tsv', sep='\t', index=False)
    pd.DataFrame(priors).to_csv(OUT / 'reference_sampling_and_assigned_fractions.tsv', sep='\t', index=False)
    pd.DataFrame(heldout_metrics).to_csv(OUT / 'heldout_gene_recovery.tsv', sep='\t', index=False)
    pd.concat(grid_values, ignore_index=True).to_csv(OUT / 'heldout_gene_grid250um_evaluation.tsv.gz', sep='\t', index=False)
    pd.DataFrame(group_markers).to_csv(OUT / 'fine_group_measured_marker_coverage.tsv', sep='\t', index=False)
    print(pd.DataFrame(summaries).to_string(index=False), flush=True)
    print(pd.DataFrame(heldout_metrics).groupby('stage')[['cell_spearman_rho', 'depth_area_gene_count_adjusted_r', 'grid250um_spearman_rho']].median().to_string(), flush=True)


if __name__ == '__main__':
    main()
