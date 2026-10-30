"""Independent ID, coordinate, sampling, prediction and filtering checks."""
from collections import Counter
import json

import anndata as ad
import numpy as np
import pandas as pd
from scipy.stats import pearsonr

from run_cytospace_xenium import SOURCE, OUT, STAGES, RUNS, digest


def main():
    meta = pd.read_csv(SOURCE / 'reference_metadata.tsv.gz', sep='\t', index_col=0)
    panel = np.load(SOURCE / 'reference_shared_panel_counts.npz')
    genes = panel['genes'].astype(str)
    # NpzFile lookups decompress a whole array: materialize each layer once.
    reference_counts = {name: panel[name] for name in ['corrected', 'original']}
    catalog = pd.read_csv(SOURCE / 'reference_group_catalog.tsv', sep='\t').set_index('microenvironment_cluster')
    recover = pd.read_csv(OUT / 'reference_panel_benchmark/panel_recoverability.tsv', sep='\t')
    markers = pd.read_csv(OUT / 'snRNA_reference_marker_values_NOT_Xenium_measurements.tsv.gz', sep='\t', index_col=0)
    validation = {'stage_checks': [], 'sampled_pair_correlations': 0, 'sampled_label_consensus_checks': 0,
                  'total_assignment_rows_checked': 0, 'reference_prediction_values_checked': 0}
    for stage in STAGES:
        folder = OUT / stage
        data = ad.read_h5ad(SOURCE / stage / f'{stage}_xenium_counts_and_candidate_mapping.h5ad')
        out = pd.read_csv(folder / 'cytospace_cell_annotations.tsv.gz', sep='\t', index_col=0, float_precision='round_trip')
        assert np.array_equal(data.obs_names, out.index)
        assert np.array_equal(data.obsm['spatial'], out[['x_centroid', 'y_centroid']].to_numpy())
        assert 'Kcnd2' not in data.var_names and 'Irx5' not in data.var_names
        expected = out.cytospace_eligible
        assert out.loc[expected, 'pass_expression_QC'].all() and out.loc[expected, 'pass_transfer_information'].all()
        assert not out.loc[expected, 'CM_EC_program_overlap_review_flag'].any()
        assert not out.loc[~expected, 'cyto_corrected_seed42_reference_cell_id'].notna().any()
        # AnnData0.13 can expose base X under the sentinel layer key None.
        assert all('predict' not in str(key).lower() for key in data.layers.keys())
        source_sha = digest(SOURCE / stage / f'{stage}_xenium_counts_and_candidate_mapping.h5ad')
        query = data[:, genes].X.toarray()
        assignments = {}
        for run in RUNS:
            dest = folder / run
            manifest = json.loads((dest / 'run_manifest.json').read_text())
            assert manifest['source_query_sha256'] == source_sha
            assert manifest['upstream_source_modified'] is False
            assert manifest['assignments_sha256'] == digest(dest / 'assigned_locations.csv.gz')
            a = pd.read_csv(dest / 'assigned_locations.csv.gz', float_precision='round_trip')
            assert a.SpotID.is_unique and set(a.SpotID) == set(out.index[expected])
            assert len(a) == manifest['query_cells']
            oi = out.index.get_indexer(a.SpotID)
            assert np.array_equal(a[['x_centroid', 'y_centroid']].to_numpy(), out.iloc[oi][['x_centroid', 'y_centroid']].to_numpy())
            assert np.array_equal(a.CellType.to_numpy(), out.iloc[oi].cytospace_input_lineage.to_numpy())
            ri = meta.index.get_indexer(a.OriginalCID)
            assert (ri >= 0).all() and meta.iloc[ri].age.eq(stage).all()
            assert np.array_equal(a.OriginalCID.to_numpy(), out.iloc[oi][f'cyto_{run}_reference_cell_id'].to_numpy())
            pool = pd.read_csv(dest / 'reference_pool.tsv.gz', sep='\t', index_col=0)
            assert set(a.OriginalCID).issubset(pool.index)
            for broad, group in a.groupby('CellType'):
                available = pool.index[pool.cytospace_input_lineage.eq(broad)]
                if len(group) >= len(available):
                    assert set(available) == set(group.OriginalCID)
                else:
                    assert group.OriginalCID.is_unique
            assert set(manifest['mapping_genes']).isdisjoint(manifest['withheld_genes'])
            assert '--single-cell' in manifest['command'] and '--downsample-off' in manifest['command']
            training = np.isin(genes, manifest['mapping_genes'])
            count_key = 'original' if run.startswith('original') else 'corrected'
            annotation = pd.read_csv(dest / 'annotated_assignments.tsv.gz', sep='\t', float_precision='round_trip').set_index('SpotID')
            chosen = np.random.default_rng(59).choice(len(a), size=min(25, len(a)), replace=False)
            for index in chosen:
                q = query[oi[index], training].astype(float)
                r = reference_counts[count_key][ri[index], training].astype(float)
                q = np.log2(1 + 1e6 * q / q.sum())
                r = np.log2(1 + 1e6 * r / r.sum())
                rho = pearsonr(q, r).statistic
                assert np.isclose(rho, annotation.loc[a.SpotID.iloc[index], 'mapping_gene_pearson_r'], atol=1e-10)
                validation['sampled_pair_correlations'] += 1
            validation['total_assignment_rows_checked'] += len(a)
            assignments[run] = a.set_index('SpotID')
        vent = out.spatial_lineage_candidate.eq('Ventricular-enriched CM')
        selected = np.random.default_rng(52).choice(out.index[vent], size=100, replace=False)
        for cell in selected:
            for kind in ['family', 'group']:
                labels = [out.loc[cell, f'cyto_corrected_seed{seed}_{kind}'] for seed in [42, 43, 44]]
                count = Counter(labels)
                expected_mode = sorted(count, key=lambda x: (-count[x], x))[0]
                assert out.loc[cell, f'cyto_seed_consensus_{kind}'] == expected_mode
                assert out.loc[cell, f'cyto_seed_consensus_{kind}_votes'] == count[expected_mode]
                validation['sampled_label_consensus_checks'] += 1
        for kind, key in [('family', 'mapping_family'), ('group', 'ventricular_group')]:
            field = f'cyto_stable_{kind}_candidate'
            passed = vent & ~out[field].eq('Not supported by all checks')
            allowed = set(recover.loc[recover.stage.eq(stage) & recover.resolution.eq(key) & recover.panel_recoverable, 'label'])
            assert set(out.loc[passed, field]).issubset(allowed)
            assert (out.loc[passed, f'cyto_seed_consensus_{kind}_votes'] >= 2).all()
            assert out.loc[passed, field].eq(out.loc[passed, f'cyto_corrected_seed42_{kind}']).all()
            for alternative in ['original_seed42', 'library1_seed42', 'library2_seed42']:
                assert out.loc[passed, field].eq(out.loc[passed, f'cyto_{alternative}_{kind}']).all()
            assert (out.loc[passed, 'cyto_corrected_seed42_pearson_r'] >= .30).all()
            if kind == 'group':
                assert (out.loc[passed, 'shared_panel_transcripts'] >= 100).all()
                assert (out.loc[passed, 'cyto_consensus_group_marker_genes_detected'] >= 3).all()
                assert out.loc[passed, field].map(catalog.mapping_family).eq(out.loc[passed, 'cyto_stable_family_candidate']).all()
        for gene in ['Kcnd2', 'Irx5', 'Irx3', 'Etv1', 'Mpped2', 'Kcnip2']:
            values = []
            for seed in [42, 43, 44]:
                ids = assignments[f'corrected_seed{seed}'].loc[out.index[vent], 'OriginalCID']
                values.append(markers.loc[ids, f'{gene}_corrected_log1pCP10k'].to_numpy())
            values = np.array(values).T
            assert np.allclose(values.mean(1), out.loc[vent, f'predicted_reference_{gene}_mean_log1pCP10k'])
            assert np.allclose(np.ptp(values, axis=1), out.loc[vent, f'predicted_reference_{gene}_range_across_3seeds'])
            validation['reference_prediction_values_checked'] += values.size
        validation['stage_checks'].append({'stage': stage, 'source_cells_preserved': len(out), 'all_7_runs_complete': True,
            'coordinates_exactly_conserved': True, 'Kcnd2_Irx5_unassayed': True})
        del data, query
    # Every technical benchmark must use disjoint reference libraries.
    for stage in STAGES:
        for fold in [1, 2]:
            folder = OUT / 'reference_panel_benchmark' / stage / f'query_library{fold}'
            a = pd.read_csv(folder / 'assigned_locations.csv.gz')
            m = json.loads((folder / 'benchmark_manifest.json').read_text())
            assert a.SpotID.is_unique
            assert meta.loc[a.SpotID, 'sample'].eq(m['query_library']).all()
            assert not meta.loc[a.OriginalCID, 'sample'].eq(m['query_library']).any()
    validation['status'] = 'passed'
    (OUT / 'validation.json').write_text(json.dumps(validation, indent=2))
    print(json.dumps(validation, indent=2))


if __name__ == '__main__':
    main()
