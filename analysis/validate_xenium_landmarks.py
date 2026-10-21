"""Check matched control membership and brute-force landmark distances."""
from pathlib import Path
import json

import anndata as ad
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/xenium_ventricular_mapping'
pairs = pd.read_csv(OUT / 'landmark_matched_control_pairs.tsv.gz', sep='\t')
summary = pd.read_csv(OUT / 'landmark_depth_area_matched_summary.tsv', sep='\t')
fields = {'Slit2-associated': 'measured_Slit2_program_positive', 'Conduction overlap': 'measured_conduction_overlap_positive', 'Ankrd1-associated': 'measured_Ankrd1_program_positive'}
rng = np.random.default_rng(42)
checks = []
for stage in ['P0', 'P7', 'P14', 'P21']:
    data = ad.read_h5ad(OUT / stage / f'{stage}_xenium_counts_and_candidate_mapping.h5ad', backed='r')
    obs = data.obs
    rows = pairs[pairs.stage.eq(stage)]
    for program, group in rows.groupby('program'):
        p = obs.loc[group.positive_cell_id]
        c = obs.loc[group.control_cell_id]
        field = fields[program]
        assert p[field].all() and not c[field].any()
        assert p.spatial_lineage_candidate.eq('Ventricular-enriched CM').all() and c.spatial_lineage_candidate.eq('Ventricular-enriched CM').all()
        assert np.array_equal(p.missing_nucleus_area_flag, c.missing_nucleus_area_flag)
        assert group.groupby('positive_cell_id').size().eq(5).all()
        expected = summary[summary.stage.eq(stage) & summary.program.eq(program)].iloc[0]
        assert abs(float(c.distance_to_Npr3_endothelial_landmark_um.median()) - expected.matched_control_median_endo_distance_um) < .001
        assert expected.unique_matched_controls == group.control_cell_id.nunique()
        assert expected.positive_segments == group.positive_cell_id.nunique()
    selected = rng.choice(len(obs), 50, replace=False)
    query = obs.iloc[selected][['x_centroid', 'y_centroid']].to_numpy()
    for field, distance in [('candidate_Npr3_endothelial_landmark', 'distance_to_Npr3_endothelial_landmark_um'), ('candidate_Upk3b_Wt1_surface_landmark', 'distance_to_Upk3b_Wt1_surface_landmark_um')]:
        anchors = obs.loc[obs[field], ['x_centroid', 'y_centroid']].to_numpy()
        truth = cdist(query, anchors).min(axis=1)
        assert np.allclose(truth, obs.iloc[selected][distance], atol=1e-8)
    checks.append({'stage': stage, 'matched_control_draws_checked': len(rows), 'brute_force_distance_cells': len(selected), 'landmark_sets': 2})
    data.file.close()
assert summary.max_abs_standardized_covariate_mean_difference.max() < .013
(OUT / 'landmark_validation.json').write_text(json.dumps({'status': 'passed', 'checks': checks, 'total_matched_control_draws': len(pairs), 'notes': 'Control draws are with replacement and are not biological replicates. These are numerical checks, not independent spatial/functional verification.'}, indent=2))
print(f'PASS: {len(pairs):,} matched controls have correct program/section/nucleus-status membership; 200 cells ×2 landmark distances checked by brute force.')
