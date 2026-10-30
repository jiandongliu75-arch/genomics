"""Describe retrieved author labels and compare an independent default run."""
import argparse
import json
import time

import anndata as ad
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from acquire_author_heart_atlas import OUT, PREVIOUS, STAGES
from run_author_parameter_baseline import VENT_STATES

CM_STATES = VENT_STATES + ['a.CM', 'av.CM']
MARKERS = 'Tnnt2 Myl2 Myl4 Myl7 Tbx3 Slit2 Sema3a Gja5 Cntn2 Ankrd1 Tlr4 Ddc Hmgcs2 Mki67 Top2a Cdca2 Lrtm1 Ppara'.split()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--wait-for-runs', action='store_true')
    args = parser.parse_args()
    stage_summary, profiles, distances, reviews, comparisons = [], [], [], [], []
    for stage in STAGES:
        folder = OUT / stage
        obs = pd.read_csv(folder / 'original_cells_with_author_labels.tsv.gz', sep='\t', index_col=0, float_precision='round_trip')
        data = ad.read_h5ad(PREVIOUS / stage / f'{stage}_xenium_counts_and_candidate_mapping.h5ad')
        assert np.array_equal(obs.index, data.obs_names)
        present = obs.present_in_author_atlas
        cm = obs.author_atlas_cell_state.isin(CM_STATES)
        ventricular = obs.author_atlas_cell_state.isin(VENT_STATES)
        record = {'stage': stage, 'original_segments': len(obs), 'published_atlas_segments': int(present.sum()),
            'not_in_published_atlas': int((~present).sum()), 'author_CM_segments': int(cm.sum()),
            'author_ventricular_segments': int(ventricular.sum()),
            'previous_ventricular_candidates': int(obs.spatial_lineage_candidate.eq('Ventricular-enriched CM').sum()),
            'shared_previous_and_author_ventricular': int((ventricular & obs.spatial_lineage_candidate.eq('Ventricular-enriched CM')).sum()),
            'author_ventricular_with_previous_CM_EC_review': int((ventricular & obs.CM_EC_program_overlap_review_flag).sum())}
        stage_summary.append(record)
        pd.crosstab(obs.author_atlas_cell_state, obs.spatial_lineage_candidate).to_csv(folder / 'author_state_vs_previous_lineage.tsv', sep='\t')
        review = obs.CM_EC_program_overlap_review_flag
        for label, n in obs.loc[review, 'author_atlas_cell_state'].value_counts().items():
            reviews.append({'stage': stage, 'previous_CM_EC_review_cells': int(review.sum()), 'author_state': label, 'segments': int(n)})
        totals = np.asarray(data.X.sum(axis=1)).ravel()
        genes = [g for g in MARKERS if g in data.var_names]
        values = data[:, genes].X.toarray()
        for label in CM_STATES:
            select = obs.author_atlas_cell_state.eq(label).to_numpy()
            if not select.any():
                continue
            for j, gene in enumerate(genes):
                counts = values[select, j]
                normal = np.log1p(np.divide(10000*counts, totals[select], out=np.zeros(len(counts)), where=totals[select] > 0))
                profiles.append({'stage': stage, 'author_CM_state': label, 'segments': int(select.sum()), 'gene': gene,
                    'measured_positive_fraction': float((counts > 0).mean()), 'mean_measured_count': float(counts.mean()),
                    'median_measured_count': float(np.median(counts)), 'mean_measured_log1pCP10k': float(normal.mean())})
        # Descriptive spatial context, not independent validation or wall depth.
        coords = obs[['x_centroid', 'y_centroid']].to_numpy()
        for landmark in ['endocardial.EC', 'arteriole.EC', 'cap.EC']:
            select = obs.author_atlas_cell_state.eq(landmark).to_numpy()
            if not select.any():
                continue
            d, _ = cKDTree(coords[select]).query(coords[ventricular], k=1)
            labels = obs.loc[ventricular, 'author_atlas_cell_state'].to_numpy()
            for state in VENT_STATES:
                selected = labels == state
                if selected.any():
                    distances.append({'stage': stage, 'author_CM_state': state, 'landmark_author_state': landmark,
                        'CM_segments': int(selected.sum()), 'landmark_segments': int(select.sum()),
                        'median_nearest_distance_um': float(np.median(d[selected])), 'all_ventricular_median_um': float(np.median(d))})
        baseline = folder / 'independent_default_parameters'
        if args.wait_for_runs:
            start = time.monotonic()
            if not (baseline / 'run_manifest.json').exists():
                print('WAIT_PARAMETER_BASELINE', stage, flush=True)
            while not (baseline / 'run_manifest.json').exists():
                if time.monotonic() - start > 1800:
                    raise TimeoutError(stage)
                time.sleep(5)
        pairs = pd.read_csv(baseline / 'annotated_assignments.tsv.gz', sep='\t', float_precision='round_trip').set_index('SpotID')
        assert pairs.index.is_unique and pairs.index.isin(obs.index).all()
        pairs['author_coarse_program'] = pairs.published_author_state.map({'v.CM': 'Working/metabolic', 'Ddc.CM': 'Working/metabolic',
            'Ankrd1.CM': 'Ankrd1-associated', 'Slit2.CM': 'Slit2-associated', 'p.CM': 'Cycling-associated'})
        pairs['reference_comparison_program'] = pairs.reference_mapping_family.replace({'Other working CM': 'Working/metabolic',
            'Conduction-program candidate': 'Slit2-associated'})
        # This is an explicit partial correspondence of two different ontologies,
        # not a ground-truth subtype accuracy estimate.
        pd.crosstab(pairs.author_coarse_program, pairs.reference_comparison_program).to_csv(folder / 'published_vs_independent_program_crosstab.tsv', sep='\t')
        previous = pd.read_csv(PREVIOUS.parents[0] / 'cytospace_xenium' / stage / 'original_seed42/annotated_assignments.tsv.gz', sep='\t').set_index('SpotID')
        previous = previous[previous.CellType.eq('Ventricular_CM')]
        shared = pairs.index.intersection(previous.index)
        comparisons.append({'stage': stage, 'default_parameter_matches': len(pairs),
            'median_original_reference_profile_r': float(pairs.original_reference_profile_pearson_r.median()),
            'median_within_ventricular_shuffled_r': float(pairs.within_ventricular_shuffled_profile_r.median()),
            'fraction_assigned_reference_gt1500_shared_molecules': float(pairs.reference_shared_molecules_gt1500.mean()),
            'partial_ontology_agreement_fraction': float(pairs.author_coarse_program.eq(pairs.reference_comparison_program).mean()),
            'overlap_with_previous_raw_run': len(shared),
            'native_family_agreement_on_overlapping_cells': float(pairs.loc[shared, 'reference_mapping_family'].eq(previous.loc[shared, 'reference_mapping_family']).mean()),
            'comparison_note': 'Published final labels are not independent ground truth; reference families and author states have a partial correspondence only. Cohort/seed/counts/partition effects are not isolated.'})
        print('SUMMARIZED', stage, record['author_ventricular_segments'], flush=True)
        del data, values
    pd.DataFrame(stage_summary).to_csv(OUT / 'stage_summary.tsv', sep='\t', index=False)
    pd.DataFrame(profiles).to_csv(OUT / 'measured_marker_profiles_by_author_CM_state.tsv', sep='\t', index=False)
    pd.DataFrame(distances).to_csv(OUT / 'descriptive_distances_to_author_EC_states.tsv', sep='\t', index=False)
    pd.DataFrame(reviews).to_csv(OUT / 'previous_CM_EC_review_vs_author_state.tsv', sep='\t', index=False)
    pd.DataFrame(comparisons).to_csv(OUT / 'independent_parameter_baseline_summary.tsv', sep='\t', index=False)
    counts = pd.read_csv(OUT / 'author_cell_state_counts.tsv', sep='\t')
    counts = counts[counts.author_cell_state.isin(CM_STATES)].pivot(index='author_cell_state', columns='stage', values='segments')
    counts.reindex(index=CM_STATES, columns=STAGES).fillna(0).astype(int).to_csv(OUT / 'published_CM_state_counts_by_stage.tsv', sep='\t')
    print(pd.DataFrame(comparisons).to_string(index=False), flush=True)


if __name__ == '__main__':
    main()
