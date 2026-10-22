"""Independently check barcode joins, native pairs and measured summaries."""
import json
import hashlib

import anndata as ad
import numpy as np
import pandas as pd
from scipy.stats import pearsonr

from acquire_author_heart_atlas import OUT, PREVIOUS, STAGES
from run_author_parameter_baseline import VENT_STATES


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as handle:
        for block in iter(lambda: handle.read(8*1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def main():
    metadata = pd.read_csv(PREVIOUS / 'reference_metadata.tsv.gz', sep='\t', index_col=0)
    with np.load(PREVIOUS / 'reference_shared_panel_counts.npz') as bundle:
        reference, genes, ids = bundle['original'], bundle['genes'].astype(str), bundle['cell_ids'].astype(str)
    np.testing.assert_array_equal(ids, metadata.index)
    source_reference_sha = sha256(PREVIOUS / 'reference_shared_panel_counts.npz')
    profiles = pd.read_csv(OUT / 'measured_marker_profiles_by_author_CM_state.tsv', sep='\t')
    summary = pd.read_csv(OUT / 'stage_summary.tsv', sep='\t').set_index('stage')
    records = []
    for stage in STAGES:
        folder = OUT / stage
        author = pd.read_csv(folder / 'author_atlas_annotations.tsv.gz', sep='\t', float_precision='round_trip').set_index('cell_id')
        joined = pd.read_csv(folder / 'original_cells_with_author_labels.tsv.gz', sep='\t', index_col=0, float_precision='round_trip')
        path = PREVIOUS / stage / f'{stage}_xenium_counts_and_candidate_mapping.h5ad'
        data = ad.read_h5ad(path)
        assert author.index.is_unique and author.index.isin(data.obs_names).all()
        np.testing.assert_array_equal(joined.index, data.obs_names)
        for name in ['x_centroid', 'y_centroid', 'transcript_counts', 'total_counts', 'control_probe_counts', 'control_codeword_counts']:
            np.testing.assert_array_equal(author[name], data.obs.loc[author.index, name])
            np.testing.assert_array_equal(joined[name], data.obs[name])
        np.testing.assert_array_equal(joined.loc[author.index, 'author_atlas_cell_state'], author['cytoSPACE.final.anno'])
        assert joined.present_in_author_atlas.sum() == len(author) == summary.loc[stage, 'published_atlas_segments']
        assert joined.loc[~joined.present_in_author_atlas, 'author_atlas_cell_state'].eq('Not present in published atlas').all()
        features = pd.read_csv(folder / 'atlas_features.tsv', sep='\t').atlas_feature
        assert not features.isin(['Kcnd2', 'Irx5']).any() and features.eq('TF_Irx5').any()
        assert len(set(features) & set(data.var_names)) == 403
        assert not data.var_names.isin(['Kcnd2', 'Irx5']).any()
        native = folder / 'independent_default_parameters'
        pairs = pd.read_csv(native / 'annotated_assignments.tsv.gz', sep='\t', float_precision='round_trip')
        original_pairs = pd.read_csv(native / 'assigned_locations.csv.gz', float_precision='round_trip')
        pd.testing.assert_frame_equal(pairs[original_pairs.columns], original_pairs)
        assert pairs.SpotID.is_unique and pairs.UniqueCID.is_unique
        assert pairs.CellType.eq('Ventricular_CM').all()
        query = data[:, genes].X.toarray()
        cohort = joined.author_atlas_cell_state.isin(VENT_STATES) & (query.sum(1) > 0)
        assert set(pairs.SpotID) == set(joined.index[cohort])
        qi = data.obs_names.get_indexer(pairs.SpotID)
        ri = metadata.index.get_indexer(pairs.OriginalCID)
        assert (ri >= 0).all()
        assert metadata.iloc[ri].age.eq(stage).all()
        assert metadata.iloc[ri].broad_class.eq('Ventricular-enriched CM').all()
        for key in ['sample', 'ventricular_group', 'mapping_family']:
            np.testing.assert_array_equal(pairs['reference_'+key], metadata.iloc[ri][key])
        np.testing.assert_array_equal(pairs.published_author_state, joined.iloc[qi].author_atlas_cell_state)
        np.testing.assert_array_equal(pairs[['x_centroid', 'y_centroid']], data.obs.iloc[qi][['x_centroid', 'y_centroid']])
        np.testing.assert_array_equal(pairs.reference_shared_molecules_gt1500, reference[ri].sum(1) > 1500)
        pool = pd.read_csv(native / 'reference_pool.tsv.gz', sep='\t', index_col=0)
        expected_pool = metadata.index[metadata.age.eq(stage) & metadata.broad_class.eq('Ventricular-enriched CM') & (reference.sum(1) > 0)]
        assert pool.index.is_unique and set(pool.index) == set(expected_pool)
        # Upstream includes each available nucleus once, then samples the
        # remaining rows with replacement; individual reuse is not bounded by
        # the mean ratio of query to reference cells.
        reuse = pairs.OriginalCID.value_counts().reindex(pool.index, fill_value=0)
        assert (reuse >= 1).all() and reuse.sum() == len(pairs)
        sampled = np.random.default_rng(71).choice(len(pairs), 50, replace=False)
        shuffled_ri = np.random.default_rng(149).permutation(ri)
        errors = []
        for i in sampled:
            q = query[qi[i]].astype(float)
            q = np.log2(1 + 1e6*q/q.sum())
            for r_index, column in [(ri[i], 'original_reference_profile_pearson_r'),
                                    (shuffled_ri[i], 'within_ventricular_shuffled_profile_r')]:
                r = reference[r_index].astype(float)
                r = np.log2(1 + 1e6*r/r.sum())
                errors.append(abs(pearsonr(q, r).statistic - pairs.iloc[i][column]))
        assert max(errors) < 1e-12
        total = np.asarray(data.X.sum(1)).ravel()
        checked = 0
        for row in profiles[profiles.stage.eq(stage)].itertuples():
            mask = joined.author_atlas_cell_state.eq(row.author_CM_state).to_numpy()
            counts = data[mask, [row.gene]].X.toarray().ravel()
            assert row.segments == mask.sum()
            normalized = np.log1p(np.divide(10000*counts, total[mask], out=np.zeros(len(counts)), where=total[mask] > 0))
            for actual, expected in [(np.mean(counts > 0), row.measured_positive_fraction),
                                     (counts.mean(), row.mean_measured_count),
                                     (np.median(counts), row.median_measured_count),
                                     (normalized.mean(), row.mean_measured_log1pCP10k)]:
                assert abs(actual-expected) < 1e-10
            checked += 1
        manifest = json.loads((native / 'run_manifest.json').read_text())
        assert manifest['source_query_sha256'] == sha256(path)
        assert manifest['source_reference_sha256'] == source_reference_sha
        assert manifest['assignments_sha256'] == sha256(native / 'assigned_locations.csv.gz')
        assert manifest['query_segments'] == len(pairs)
        records.append({'stage': stage, 'original_cells': len(joined), 'published_atlas_cells': len(author),
            'native_pairs_checked': len(pairs), 'reference_nuclei': len(pool),
            'sampled_correlations_checked': len(errors), 'maximum_sampled_correlation_error': max(errors),
            'marker_state_summaries_recomputed': checked, 'centroid_and_count_delta': 0,
            'source_hashes_unchanged': True})
        print('VALIDATED', stage, records[-1], flush=True)
        del data, query
    (OUT / 'validation.json').write_text(json.dumps({'status': 'PASS', 'stages': records,
        'scope': 'Data correspondence, native constraints, diagnostic calculations and measured summaries; not independent biological validation.'}, indent=2)+'\n')


if __name__ == '__main__':
    main()
