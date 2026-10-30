"""Consolidate measured evidence, review flags, landmarks and primary outputs."""
from pathlib import Path
import json
import shutil
import importlib.metadata

import anndata as ad
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/xenium_ventricular_mapping'
CACHE = Path('/workspace/heart-atlas/results/xenium_mapping_cache')
STAGES = ['P0', 'P7', 'P14', 'P21']
GENES = 'Ankrd1 Slit2 Gja5 Cntn2 Ddc Pola1 Rad51b Mki67 Top2a Tnnt2 Tnni3 Myl2 Myl4 Myl7 Thsd4 Creb5 Ncam1 Npr3 Upk3b Wt1 Rreb1 Ppara Cd36 Cpt1b'.split()
summary_rows, profile_rows, distance_rows, fine_rows = [], [], [], []
catalog = pd.read_csv(OUT / 'reference_group_catalog.tsv', sep='\t')
reference_directory = CACHE if (CACHE / 'reference_shared_panel_counts.npz').exists() else OUT
reference_metadata = pd.read_csv(reference_directory / 'reference_metadata.tsv.gz', sep='\t', index_col=0)
reference_counts = np.load(reference_directory / 'reference_shared_panel_counts.npz')

for stage in STAGES:
    folder = OUT / stage
    path = folder / f'{stage}_xenium_counts_and_candidate_mapping.h5ad'
    data = ad.read_h5ad(path)
    obs = data.obs.copy()
    assert 'spatial_mapped_family' in obs
    def counts(gene):
        return data.X[:, data.var_names.get_loc(gene)].toarray().ravel()
    # Distinguish any marker co-detection from enrichment above the section's
    # ventricular background. Neither establishes completed cardiomyocyte division.
    obs['measured_cycling_marker_codetection'] = obs.measured_cycling_genes.ge(3) & (counts('Mki67') + counts('Top2a') >= 1) & obs.spatial_lineage_candidate.eq('Ventricular-enriched CM')
    obs['measured_cycling_program_positive'] = obs.measured_cycling_marker_codetection & (counts('Mki67') + counts('Top2a') >= 2) & obs.measured_cycling_within_section_vent_z.ge(.5)
    # Diagnostic direct-transfer labels are archived separately, never presented
    # as primary cell identities after observed cross-assay contradictions.
    archive = OUT / 'unintegrated_reference_transfer_diagnostics' / stage
    archive.mkdir(parents=True, exist_ok=True)
    for file in list(folder.glob('*')):
        if file.is_file() and (file.name.startswith(('lineage_', 'family_', 'fine_group_')) or file.name in ['cell_mapping_and_QC.tsv.gz', 'mapping_counts.tsv', 'section_summary.json', 'assay_aligned_broad_diagnostic.tsv.gz']):
            shutil.move(str(file), str(archive / file.name))
    diagnostics = [column for column in obs if column.startswith(('lineage_', 'family_', 'fine_group_')) or column in ['mapped_lineage', 'mapped_family', 'mapped_fine_group', 'CM_chamber_candidate', 'CM_chamber_score']]
    diagnostics += [column for column in obs if column.startswith('assay_aligned_lineage_')]
    if diagnostics:
        obs[diagnostics].to_csv(archive / 'diagnostic_scores_and_labels.tsv.gz', sep='\t', na_rep='NA')
        obs = obs.drop(columns=diagnostics)
    coords = obs[['x_centroid', 'y_centroid']].to_numpy()
    endo = (counts('Npr3') >= 2) & obs.endothelial_regulators_detected.ge(1).to_numpy() & obs.endothelial_vascular_markers_detected.ge(2).to_numpy() & obs.pass_expression_QC.to_numpy()
    surface = (counts('Upk3b') >= 1) & (counts('Wt1') >= 1) & ~obs.measured_CM_marker_gate.to_numpy() & obs.pass_expression_QC.to_numpy()
    obs['candidate_Npr3_endothelial_landmark'] = endo
    obs['candidate_Upk3b_Wt1_surface_landmark'] = surface
    for label, mask in [('Npr3_endothelial', endo), ('Upk3b_Wt1_surface', surface)]:
        if mask.sum():
            distance = cKDTree(coords[mask]).query(coords, k=1)[0]
        else:
            distance = np.full(len(obs), np.nan)
        obs[f'distance_to_{label}_landmark_um'] = distance
    vent = obs.spatial_lineage_candidate.eq('Ventricular-enriched CM').to_numpy()
    total = np.asarray(data.X.sum(axis=1)).ravel()
    selected = data.X[:, data.var_names.get_indexer(GENES)].toarray()
    log = np.log1p(selected * 10000 / np.maximum(total[:, None], 1))
    groups = obs.loc[vent].groupby('spatial_mapped_family', observed=True).indices
    global_vent = np.flatnonzero(vent)
    for family, local in groups.items():
        ids = global_vent[local]
        for i, gene in enumerate(GENES):
            profile_rows.append({'stage': stage, 'mapped_family': family, 'gene': gene, 'segmented_cells': len(ids), 'mean_measured_count': float(selected[ids, i].mean()), 'fraction_measured_ge1': float((selected[ids, i] >= 1).mean()), 'mean_log1p_all403_normalized_10000': float(log[ids, i].mean())})
    populations = {'all_ventricular_candidates': vent}
    for field in ['measured_Ankrd1_program_positive', 'measured_Slit2_program_positive', 'measured_conduction_overlap_positive', 'measured_cycling_program_positive']:
        populations[field] = obs[field].to_numpy(dtype=bool)
    for name, mask in populations.items():
        distance_rows.append({'stage': stage, 'population': name, 'segmented_cells': int(mask.sum()), 'Npr3_landmarks': int(endo.sum()), 'Upk3b_Wt1_landmarks': int(surface.sum()), 'median_nearest_Npr3_landmark_um': float(obs.loc[mask, 'distance_to_Npr3_endothelial_landmark_um'].median()) if mask.any() else np.nan, 'median_nearest_surface_landmark_um': float(obs.loc[mask, 'distance_to_Upk3b_Wt1_surface_landmark_um'].median()) if mask.any() else np.nan, 'median_transcripts': float(obs.loc[mask, 'transcript_counts'].median()) if mask.any() else np.nan, 'median_cell_area_um2': float(obs.loc[mask, 'cell_area'].median()) if mask.any() else np.nan})
    family_assigned = vent & obs.spatial_mapped_family.ne('Unresolved ventricular state').to_numpy()
    fine_assigned = obs.spatial_mapped_fine_group.ne('Not resolved').to_numpy()
    info = {'stage': stage, 'source_segmented_cells': len(obs), 'QC_pass': int(obs.pass_expression_QC.sum()), 'measured_CM_marker_gate': int(obs.measured_CM_marker_gate.sum()), 'ventricular_CM_candidates': int(vent.sum()), 'CM_EC_overlap_review': int(obs.CM_EC_program_overlap_review_flag.sum()), 'family_assigned': int(family_assigned.sum()), 'fine_group_assigned': int(fine_assigned.sum()), 'measured_Ankrd1_program_positive': int(obs.measured_Ankrd1_program_positive.sum()), 'measured_Slit2_program_positive': int(obs.measured_Slit2_program_positive.sum()), 'measured_conduction_overlap_positive': int(obs.measured_conduction_overlap_positive.sum()), 'measured_cycling_program_positive': int(obs.measured_cycling_program_positive.sum()), 'candidate_Npr3_landmarks': int(endo.sum()), 'candidate_surface_landmarks': int(surface.sum())}
    summary_rows.append(info)
    (folder / 'spatial_summary.json').write_text(json.dumps(info, indent=2))
    check = pd.read_csv(folder / 'assay_aligned_fine_group_recoverability.tsv', sep='\t').set_index('reference_class')
    assigned_counts = obs.spatial_mapped_fine_group.value_counts()
    for row in catalog[catalog.stage.eq(stage)].itertuples():
        group = row.microenvironment_cluster
        fine_rows.append({'stage': stage, 'reference_group': group, 'reference_identity': row.candidate_identity, 'reference_support_status': row.support_status, 'reference_nuclei': row.nuclei, 'aligned_reference_panel_supported': bool(check.loc[group, 'panel_supported']), 'reference_min_heldout_precision': check.loc[group, 'min_precision'], 'reference_min_heldout_recall': check.loc[group, 'min_recall'], 'mapped_spatial_candidate_segments': int(assigned_counts.get(group, 0)), 'interpretation': 'Reference-like candidate only; zero means not resolved, not biological absence'})
    obs.to_csv(folder / 'spatial_candidates_and_measured_programs.tsv.gz', sep='\t', na_rep='NA', float_format='%.7g')
    data.obs = obs
    stage_mask = reference_metadata.age.eq(stage).to_numpy()
    stage_metadata = reference_metadata.loc[stage_mask]
    panel = reference_counts['corrected'][stage_mask]
    information = (panel.sum(axis=1) >= 10) & ((panel >= 1).sum(axis=1) >= 5)
    vr = information & stage_metadata.broad_class.eq('Ventricular-enriched CM').to_numpy()
    gr = information & stage_metadata.broad_class.ne('Unresolved/mixed').to_numpy()
    corrected_embedding = np.load(CACHE / f'{stage}_ventricular_marker_v2_corrected_joint_assay_embedding.npz')
    original_embedding = np.load(CACHE / f'{stage}_ventricular_marker_v2_original_joint_assay_embedding.npz')
    global_embedding = np.load(CACHE / f'{stage}_broad_corrected_joint_assay_embedding.npz')
    assert int(corrected_embedding['query_n']) == vent.sum()
    assert corrected_embedding['reference_embedding'].shape[0] == vr.sum()
    np.savez_compressed(folder / 'reference_aligned_embeddings.npz', ventricular_corrected=corrected_embedding['reference_embedding'], ventricular_original=original_embedding['reference_embedding'], ventricular_reference_cell_ids=stage_metadata.index[vr].to_numpy(dtype=str), global_corrected=global_embedding['reference_embedding'], global_reference_cell_ids=stage_metadata.index[gr].to_numpy(dtype=str))
    data.obsm['X_assay_aligned_ventricular_original'] = np.full((len(obs), original_embedding['query_embedding'].shape[1]), np.nan, dtype=np.float32)
    data.obsm['X_assay_aligned_ventricular_original'][vent] = original_embedding['query_embedding']
    data.uns['primary_mapping_columns'] = ['spatial_lineage_candidate', 'spatial_mapped_family', 'spatial_mapped_fine_group']
    data.uns['unintegrated_classifier_columns'] = 'Archived separately as diagnostics; not primary identities.'
    data.uns['landmark_limits'] = 'Nearest marker-positive cell distances in 2D micrometers, not traced endocardial/epicardial boundaries or normalized transmural depth.'
    data.uns['reference_embedding_file'] = 'reference_aligned_embeddings.npz; reference IDs link to root reference_metadata.tsv.gz'
    temporary = path.with_name(path.stem + '.partial.h5ad')
    data.write_h5ad(temporary, compression='gzip')
    temporary.replace(path)
    pd.DataFrame([{'stage': stage, 'level': field, 'label': label, 'candidate_segments': int(n)} for field in ['spatial_lineage_candidate', 'spatial_mapped_family', 'spatial_mapped_fine_group'] for label, n in obs[field].value_counts().items()]).to_csv(folder / 'spatial_candidate_counts.tsv', sep='\t', index=False)
    print('FINALIZED_SPATIAL', info, flush=True)

pd.DataFrame(summary_rows).to_csv(OUT / 'section_summary.tsv', sep='\t', index=False)
pd.DataFrame(profile_rows).to_csv(OUT / 'measured_gene_profiles_by_mapped_family.tsv', sep='\t', index=False)
pd.DataFrame(distance_rows).to_csv(OUT / 'candidate_landmark_distance_summary.tsv', sep='\t', index=False, na_rep='NA')
pd.DataFrame(fine_rows).to_csv(OUT / 'fine_group_spatial_mapping_coverage.tsv', sep='\t', index=False)
if reference_directory != OUT:
    shutil.copy2(reference_directory / 'reference_metadata.tsv.gz', OUT / 'reference_metadata.tsv.gz')
    shutil.copy2(reference_directory / 'reference_shared_panel_counts.npz', OUT / 'reference_shared_panel_counts.npz')
ref = pd.read_csv(ROOT / 'results/ventricular_additional_axes/electrical_axis_cluster_profiles.tsv', sep='\t')
ref['measurement_assay'] = 'snRNA-seq reference; NOT a Xenium measurement'
ref.to_csv(OUT / 'reference_electrical_profiles_NOT_spatial_measurements.tsv', sep='\t', index=False, na_rep='NA')
(OUT / 'measured_program_definitions.json').write_text(json.dumps({'Ankrd1_associated': ['Ankrd1', 'Creb5', 'Tlr4', 'Itga9', 'Pdlim3'], 'Slit2_associated': ['Slit2', 'Sema3a', 'Ephb1', 'Gja5', 'Thsd4'], 'conduction_overlap': ['Cntn2', 'Gja5', 'Robo1', 'Cpne5', 'Sema3a', 'Tbx3'], 'cycling': ['Mki67', 'Top2a', 'Cenpf', 'Cenpe', 'Cdca2', 'Kif11', 'Kif15', 'Smc4', 'Hjurp', 'Atad2'], 'DNA_maintenance': ['Pola1', 'Rad51b', 'Hells', 'Wdr76', 'Rrm1', 'Ezh2'], 'lipid_regulatory': ['Cd36', 'Cpt1b', 'Ppara', 'Ppargc1a', 'Ppargc1b'], 'scores': 'Mean log1p(10000 * measured count / exact-shared-panel total), z-scored against marker-supported ventricular candidates separately within each section. Partial selected panels, not the full author signatures or physical exposure measurements.', 'positive_flags': {'Ankrd1': 'Ankrd1>=2 molecules, >=3 panel genes, section-relative z>=0.5.', 'Slit2': 'Slit2>=1, >=2 panel genes, z>=0.5.', 'conduction_overlap': 'Cntn2>=1 and Gja5>=1 and >=3 conduction-overlap genes. Not verified Purkinje identity.', 'cycling': '>=3 cycling panel genes, Mki67+Top2a>=2 molecules, z>=0.5. Co-detection-only flag also retained. No completed-division inference.'}, 'landmarks': {'Npr3': 'Npr3>=2 with >=1 Erg/Fli1 and >=2 vascular markers, QC-pass.', 'surface': 'Upk3b>=1 and Wt1>=1, outside measured CM gate, QC-pass.', 'limits': '2D nearest marker-positive cell distances, not wall segmentation, transmural depth or replicated anatomical tests.'}}, indent=2))
methods = {
    'primary_columns': ['spatial_lineage_candidate', 'spatial_mapped_family', 'spatial_mapped_fine_group'],
    'shared_features': 381,
    'assay_alignment': 'Per-stage exact-shared-panel log1p(10000 * counts / shared total); gene means/SD fitted to balanced up-to-20k reference/query cells (SD floor0.5), clip standardized values to[-5,5]; PCA30 global context, PCA20 ventricular subset; Harmony over assay only, theta2/lambda1, nclust50(global)/30(ventricular), max20 Harmony/max10 kmeans iterations, seed42,4cores. Counts and physical coordinates unchanged.',
    'primary_CM_gate': '>=3 of Tnnt2/Tnni3/Myh6/Myh7/Myl2/Actc1, >=10 corresponding transcripts, >=10% of measured403-gene counts, CM-marker fraction>=1.5 times each other-lineage marker fraction; pass QC/transfer-information thresholds; withhold CM/EC program-overlap review cells.',
    'CM_EC_overlap': '>=1 of Erg/Fli1 plus>=3 of Cdh5/Pecam1/Kdr/Esam/Tek/Vwf/Flt1/Tie1 and CM-marker fraction>=10%. This is a conservative review category, not a doublet or endothelial identity claim. Source paper reports capillary ECs with muscle-associated RNA.',
    'ventricular_gate': 'Measured-CM-positive, Myl2>=2 and Myl2>=Myl4. Myl7 is NOT used as negative evidence because immature ventricular CMs express it. Label is ventricular-enriched candidate, not traced chamber anatomy.',
    'atrial_gate': 'Measured-CM-positive, Myl4>=2, Myl7>=1, Myl2<=1. Remaining CM-positive segments retain chamber-unresolved label.',
    'other_lineages': 'Global assay-aligned nearest-reference vote plus positive>=2 marker genes and marker fraction>=5% (epicardial2%), score>=0.65/margin>=0.15, supported heldout-reference class. These are secondary context candidates, not a full independently validated spatial atlas.',
    'transfer': '30-nearest reference neighbors, cosine distance in assay-aligned coordinates; inverse-distance weighted votes with distance floor0.01. Family/fine labels use separate ventricular-only corrected-reference and original-reference alignments.',
    'reference_benchmark': 'Hold out each original same-age reference library in turn within the joint embedding. Unsupervised embedding includes unlabeled reference/query cells (transductive benchmark). Require minprecision>=0.60,minrecall>=0.35,minF1>=0.50,>=30train/test nuclei per class. This is panel recoverability, not independently validated Xenium accuracy or animal-level inference.',
    'family_acceptance': 'Supported reference family, vote>=0.60/margin>=0.15,>=50shared-panel transcripts, same label from both original reference libraries and separate original-count-reference alignment; cycling/conduction candidates require>=2 respective measured genes. Source uncharacterized P7_M7/M8 not used as family classes.',
    'fine_acceptance': 'Supported reference group and>=3 upstream both-library raw/corrected markers, vote>=0.70/margin>=0.20,>=100shared-panel transcripts, same label from both libraries and original-count-reference alignment, consistent with accepted family. No forced label for other cells.',
    'QC': json.loads((OUT / 'mapping_methods.json').read_text())['query_QC'],
    'diagnostics': 'Unintegrated logistic transfer and global chamber-vote diagnostics archived separately after cross-assay contradictions. Those labels are not primary spatial identities. Low100-molecule benchmark remains a coverage sensitivity, not a claim of actual Xenium accuracy.',
    'spatial': 'Original µm centroids registered to acquired DAPI using experiment pixel_size0.2125µm at full resolution. Original segmentation polygons used in illustrative zooms. No spatial label smoothing, anatomical warping or transcript reassignment.',
    'Kcnd2': 'Not present by gene symbol OR EnsemblID ENSMUSG00000060882 in supplied panel; no measurement, imputation or Kcnd2-derived orientation performed.',
    'Irx5': 'Not present by gene symbol OR EnsemblID ENSMUSG00000031737 in supplied panel.',
    'limitations': ['One deposited acquisition region per stage analyzed; spatial biological replication is not established from these four regions.', 'Markers/thresholds and expression transfer can be affected by cell fragments, segmentation leakage, sparse detection and true cross-lineage RNA programs.', 'Partial marker programs and reference matches are hypotheses, not functional, chamber, transmural or Purkinje verification.', 'The same measured genes influence mapping and marker-consistency plots; those plots are not independent validation.', 'Gene normalization over the targeted panel is compositional; gene-specific probe efficiencies and nuclear versus cellular RNA differences remain despite alignment.', 'Counts of accepted matches are operational filtered segment counts, not subtype prevalence or CM proliferation rates.', 'Nearest-landmark distances are 2D marker-positive-cell distances, not wall depth or causal effects.'],
    'versions': {package: importlib.metadata.version(package) for package in ['numpy', 'scipy', 'pandas', 'anndata', 'scikit-learn', 'harmonypy', 'h5py', 'matplotlib', 'tifffile', 'imagecodecs', 'remotezip']}
}
(OUT / 'primary_spatial_methods.json').write_text(json.dumps(methods, indent=2))
