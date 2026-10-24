"""Descriptive landmark proximity with RNA-depth/area matched CM controls."""
from pathlib import Path
import os
import json

os.environ.setdefault('MPLCONFIGDIR', '/tmp/heart-mpl')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/xenium_ventricular_mapping'
STAGES = ['P0', 'P7', 'P14', 'P21']
PROGRAMS = [('Slit2-associated', 'measured_Slit2_program_positive'), ('Conduction overlap', 'measured_conduction_overlap_positive'), ('Ankrd1-associated', 'measured_Ankrd1_program_positive')]
rows, pairs = [], []
for stage in STAGES:
    obs = pd.read_csv(OUT / stage / 'spatial_candidates_and_measured_programs.tsv.gz', sep='\t', index_col=0)
    obs = obs[obs.spatial_lineage_candidate.eq('Ventricular-enriched CM')].copy()
    cov = np.column_stack([np.log1p(obs.transcript_counts), np.log1p(obs.cell_area), np.log1p(obs.measured_gene_count)])
    cov = (cov - cov.mean(axis=0)) / np.maximum(cov.std(axis=0), 1e-8)
    for name, field in PROGRAMS:
        positive = obs[field].to_numpy(dtype=bool)
        selected_positive, selected_control = [], []
        for nucleus_status in [False, True]:
            stratum = obs.missing_nucleus_area_flag.eq(nucleus_status).to_numpy()
            pos = np.flatnonzero(positive & stratum)
            control = np.flatnonzero(~positive & stratum)
            if not len(pos) or not len(control):
                continue
            model = NearestNeighbors(n_neighbors=min(5, len(control)), algorithm='kd_tree').fit(cov[control])
            distance, neighbor = model.kneighbors(cov[pos])
            for j, p in enumerate(pos):
                for k, c in enumerate(control[neighbor[j]]):
                    selected_positive.append(p)
                    selected_control.append(c)
                    pairs.append({'stage': stage, 'program': name, 'positive_cell_id': obs.index[p], 'control_cell_id': obs.index[c], 'control_slot': k + 1, 'standardized_nuisance_distance': float(distance[j, k])})
        p = np.array(selected_positive, dtype=int)
        c = np.array(selected_control, dtype=int)
        assert len(p) and len(p) == len(c)
        rows.append({'stage': stage, 'program': name, 'positive_segments': int(positive.sum()), 'controls_available': int((~positive).sum()), 'matched_control_draws': len(c), 'unique_matched_controls': len(np.unique(c)), 'positive_median_endo_distance_um': float(obs.loc[positive, 'distance_to_Npr3_endothelial_landmark_um'].median()), 'negative_median_endo_distance_um': float(obs.loc[~positive, 'distance_to_Npr3_endothelial_landmark_um'].median()), 'matched_control_median_endo_distance_um': float(obs.iloc[c].distance_to_Npr3_endothelial_landmark_um.median()), 'positive_median_surface_distance_um': float(obs.loc[positive, 'distance_to_Upk3b_Wt1_surface_landmark_um'].median()), 'matched_control_median_surface_distance_um': float(obs.iloc[c].distance_to_Upk3b_Wt1_surface_landmark_um.median()), 'positive_median_transcripts': float(obs.iloc[p].transcript_counts.median()), 'matched_median_transcripts': float(obs.iloc[c].transcript_counts.median()), 'positive_median_cell_area_um2': float(obs.iloc[p].cell_area.median()), 'matched_median_cell_area_um2': float(obs.iloc[c].cell_area.median()), 'max_abs_standardized_covariate_mean_difference': float(np.abs(cov[p].mean(axis=0) - cov[c].mean(axis=0)).max())})
table = pd.DataFrame(rows)
table.to_csv(OUT / 'landmark_depth_area_matched_summary.tsv', sep='\t', index=False)
pd.DataFrame(pairs).to_csv(OUT / 'landmark_matched_control_pairs.tsv.gz', sep='\t', index=False)
fig, axes = plt.subplots(1, 3, figsize=(16, 5.7), sharey=True)
for ax, (name, _) in zip(axes, PROGRAMS):
    sub = table[table.program.eq(name)].set_index('stage').loc[STAGES]
    for i, (_, row) in enumerate(sub.iterrows()):
        ax.plot([i - .10, i + .10], [row.positive_median_endo_distance_um, row.matched_control_median_endo_distance_um], color='#999999', lw=1.5)
    ax.scatter(np.arange(4) - .10, sub.positive_median_endo_distance_um, color='#8e44ad', label='Program-positive', s=55)
    ax.scatter(np.arange(4) + .10, sub.matched_control_median_endo_distance_um, color='#666666', marker='s', label='Matched negative controls', s=45)
    ax.set_xticks(range(4), STAGES)
    ax.set_title(name)
    ax.set_ylim(bottom=0)
    ax.spines[['top', 'right']].set_visible(False)
axes[0].set_ylabel('Median nearest Npr3-positive endothelial landmark (µm)')
axes[-1].legend(frameon=False, fontsize=9, loc='upper left')
fig.suptitle('Descriptive spatial association after matching RNA depth and cell area', fontsize=16)
fig.text(.5, .025, 'Within-section controls matched on log transcript count, cell area and detected genes; missing-nucleus status matched exactly.\nFive nearest controls with replacement per positive segment; no spatial coordinates used for matching.\nPoints summarize one deposited region per age, without inferential error bars or animal-level tests. Landmarks are not traced wall depth.', ha='center', fontsize=9, linespacing=1.45)
fig.subplots_adjust(left=.065, right=.98, bottom=.25, top=.85, wspace=.15)
fig.savefig(OUT / 'landmark_depth_area_matched.png', dpi=190, facecolor='white')
fig.savefig(OUT / 'landmark_depth_area_matched.pdf', facecolor='white')
plt.close(fig)
(OUT / 'landmark_matching_methods.json').write_text(json.dumps({'comparison': 'Measured program-positive versus program-negative marker-supported ventricular segments within the SAME section.', 'matching': 'Five nearest negative controls with replacement in standardized log1p measured transcript count, log1p cell area and log1p detected genes; exact strata for missing-nucleus-area status.', 'excluded_from_matching': 'No spatial coordinates, landmark distances or subtype labels used as matching covariates.', 'limits': 'Descriptive conditional sensitivity, not causal inference or an independent-animal test. RNA depth and cell area can themselves be biological mediators. Marker landmarks and segmentation need anatomical validation. Control draws are not independent replicates.'}, indent=2))
print(table[['stage', 'program', 'positive_median_endo_distance_um', 'matched_control_median_endo_distance_um', 'max_abs_standardized_covariate_mean_difference']].round(3).to_string(index=False))
