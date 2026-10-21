"""Assay-adjusted spatial candidates with measured CM gates and library checks.

Harmony changes coordinates only. All published Xenium expression remains original
integer counts. The unintegrated logistic analysis is a diagnostic comparator.
"""
from pathlib import Path
import argparse
import json
import time

import anndata as ad
import harmonypy
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn import config_context
from sklearn.decomposition import PCA
from sklearn.metrics import classification_report
from sklearn.neighbors import NearestNeighbors

ROOT = Path(__file__).resolve().parents[1]
BASE = Path('/workspace/heart-atlas')
OUT = ROOT / 'results/xenium_ventricular_mapping'
CACHE = BASE / 'results/xenium_mapping_cache'
MARKERS = {
    'CM': 'Tnnt2 Tnni3 Myh6 Myh7 Myl2 Actc1'.split(),
    'Endothelial': 'Pecam1 Cdh5 Kdr Esam Tek Vwf Flt1 Tie1'.split(),
    'Fibroblast': 'Col1a1 Col1a2 Col3a1 Dcn Lum Pdgfra Tcf21'.split(),
    'Pericyte': 'Rgs5 Cspg4 Pdgfrb'.split(),
    'Smooth muscle': 'Myh11 Acta2 Tagln Cnn1'.split(),
    'Immune': 'Ptprc Csf1r Mertk Mrc1 Fyb Cd22 Cd28 Cd8b1'.split(),
    'Epicardial': 'Upk3b Wt1'.split(),
}
PROGRAMS = {
    'Ankrd1_associated': ['Ankrd1', 'Creb5', 'Tlr4', 'Itga9', 'Pdlim3'],
    'Slit2_associated': ['Slit2', 'Sema3a', 'Ephb1', 'Gja5', 'Thsd4'],
    'conduction_overlap': ['Cntn2', 'Gja5', 'Robo1', 'Cpne5', 'Sema3a', 'Tbx3'],
    'cycling': ['Mki67', 'Top2a', 'Cenpf', 'Cenpe', 'Cdca2', 'Kif11', 'Kif15', 'Smc4', 'Hjurp', 'Atad2'],
    'DNA_maintenance': ['Pola1', 'Rad51b', 'Hells', 'Wdr76', 'Rrm1', 'Ezh2'],
    'lipid_regulatory': ['Cd36', 'Cpt1b', 'Ppara', 'Ppargc1a', 'Ppargc1b'],
}


def normalize(counts):
    counts = np.asarray(counts, dtype=np.float32)
    return np.log1p(counts * (10000 / np.maximum(counts.sum(axis=1), 1))[:, None])


def embedding(reference, query, key, components=30, clusters=50):
    path = CACHE / f'{key}_joint_assay_embedding.npz'
    if path.exists():
        saved = np.load(path)
        assert int(saved['reference_n']) == len(reference) and int(saved['query_n']) == len(query)
        return saved['reference_embedding'], saved['query_embedding']
    expression = np.concatenate([normalize(reference), normalize(query)], axis=0)
    nref = len(reference)
    rng = np.random.default_rng(42)
    n = min(20000, nref, len(query))
    balanced = np.concatenate([rng.choice(nref, n, replace=False), nref + rng.choice(len(query), n, replace=False)])
    mean = expression[balanced].mean(axis=0)
    scale = np.maximum(expression[balanced].std(axis=0), .5)
    expression = np.clip((expression - mean) / scale, -5, 5)
    pca = PCA(n_components=components, svd_solver='randomized', random_state=42)
    pca.fit(expression[balanced])
    pc = pca.transform(expression).astype(np.float32)
    del expression
    meta = pd.DataFrame({'assay': ['snRNA'] * nref + ['Xenium'] * len(query)})
    print(f'ALIGN {key} {pc.shape}', flush=True)
    fitted = harmonypy.run_harmony(pc, meta, ['assay'], theta=2, lamb=1, nclust=clusters, max_iter_harmony=20, max_iter_kmeans=10, random_state=42, ncores=4, verbose=False)
    aligned = np.asarray(fitted.Z_corr)
    aligned = aligned if aligned.shape == pc.shape else aligned.T
    assert aligned.shape == pc.shape and np.isfinite(aligned).all()
    np.savez_compressed(path, reference_embedding=aligned[:nref], query_embedding=aligned[nref:], reference_n=nref, query_n=len(query), pca=pc, gene_mean=mean, gene_scale=scale, pca_components=pca.components_)
    return aligned[:nref], aligned[nref:]


def vote(reference_embedding, query_embedding, labels, classes=None, k=30):
    labels = np.asarray(labels, dtype=str)
    classes = np.unique(labels) if classes is None else np.asarray(classes, dtype=str)
    codes = pd.Index(classes).get_indexer(labels)
    assert (codes >= 0).all()
    model = NearestNeighbors(n_neighbors=min(k, len(labels)), metric='cosine', algorithm='brute', n_jobs=4)
    model.fit(reference_embedding)
    with config_context(working_memory=128):
        distances, neighbors = model.kneighbors(query_embedding)
    weights = 1 / np.maximum(distances, .01)
    scores = np.zeros((len(query_embedding), len(classes)), dtype=np.float32)
    for i in range(neighbors.shape[1]):
        np.add.at(scores, (np.arange(len(scores)), codes[neighbors[:, i]]), weights[:, i])
    scores /= scores.sum(axis=1)[:, None]
    order = np.argsort(scores, axis=1)
    best = classes[order[:, -1]]
    top = scores[np.arange(len(scores)), order[:, -1]]
    margin = top - scores[np.arange(len(scores)), order[:, -2]]
    return best, top, margin, distances[:, 0], scores, classes


def benchmark(stage, level, representation, labels, libraries, folder):
    labels = np.asarray(labels, dtype=str)
    libraries = np.asarray(libraries, dtype=str)
    classes = np.unique(labels)
    rows = []
    for library in sorted(np.unique(libraries)):
        test = libraries == library
        best, _, _, _, _, _ = vote(representation[~test], representation[test], labels[~test], classes)
        report = classification_report(labels[test], best, labels=classes, output_dict=True, zero_division=0)
        for label in classes:
            item = report[label]
            rows.append({'stage': stage, 'level': level, 'heldout_library': library, 'reference_class': label, 'precision': item['precision'], 'recall': item['recall'], 'f1': item['f1-score'], 'test_nuclei': int(item['support']), 'training_nuclei': int((labels[~test] == label).sum())})
    table = pd.DataFrame(rows)
    table.to_csv(folder / f'assay_aligned_{level}_heldout_metrics.tsv', sep='\t', index=False)
    summary = table.groupby('reference_class').agg(min_precision=('precision', 'min'), min_recall=('recall', 'min'), min_f1=('f1', 'min'), min_test_nuclei=('test_nuclei', 'min'), min_training_nuclei=('training_nuclei', 'min'))
    summary['panel_supported'] = summary.min_precision.ge(.60) & summary.min_recall.ge(.35) & summary.min_f1.ge(.50) & summary.min_test_nuclei.ge(30) & summary.min_training_nuclei.ge(30)
    summary.to_csv(folder / f'assay_aligned_{level}_recoverability.tsv', sep='\t')
    print(f'{stage} ALIGNED {level} SUPPORTED {summary.index[summary.panel_supported].tolist()}', flush=True)
    return summary


def fraction_and_detection(matrix, measured, genes):
    ids = [measured.index(gene) for gene in genes if gene in measured]
    counts = np.asarray(matrix[:, ids].sum(axis=1)).ravel()
    detected = np.asarray((matrix[:, ids] >= 1).sum(axis=1)).ravel()
    total = np.asarray(matrix.sum(axis=1)).ravel()
    return counts / np.maximum(total, 1), detected, counts


def process(stage, cache, metadata, sharedgenes, only_broad=False, reuse_global=False):
    start = time.monotonic()
    folder = OUT / stage
    data = ad.read_h5ad(folder / f'{stage}_xenium_counts_and_candidate_mapping.h5ad')
    obs = data.obs.copy()
    previous_lineage = obs.spatial_lineage_candidate.astype(str).copy() if reuse_global else None
    previous_CM_gate = obs.measured_CM_marker_gate.to_numpy().copy() if reuse_global else None
    measured = list(data.var_names)
    query = data.X[:, [measured.index(gene) for gene in sharedgenes]].toarray().astype(np.float32)
    stage_ref = metadata.age.eq(stage).to_numpy()
    corrected, original = cache['corrected'][stage_ref], cache['original'][stage_ref]
    ref = metadata.loc[stage_ref].copy()
    information = (corrected.sum(axis=1) >= 10) & ((corrected >= 1).sum(axis=1) >= 5)
    keepref = information & ref.broad_class.ne('Unresolved/mixed').to_numpy()
    keepquery = obs.pass_transfer_information.to_numpy()
    ri = np.flatnonzero(keepref)
    qi = np.flatnonzero(keepquery)
    re, qe = embedding(corrected[ri], query[qi], f'{stage}_broad_corrected', 30, 50)
    broadcheck = benchmark(stage, 'lineage', re, ref.iloc[ri].broad_class, ref.iloc[ri]['sample'], folder)
    if reuse_global:
        best = obs.iloc[qi].assay_aligned_lineage_candidate.to_numpy(dtype=str)
        score = obs.iloc[qi].assay_aligned_lineage_score.to_numpy()
        margin = np.full(len(qi), np.nan)
    else:
        best, score, margin, _, _, _ = vote(re, qe, ref.iloc[ri].broad_class)
    obs['assay_aligned_lineage_candidate'] = pd.NA
    obs['assay_aligned_lineage_score'] = np.nan
    obs['assay_aligned_lineage_candidate'] = obs.assay_aligned_lineage_candidate.astype(object)
    obs.loc[obs.index[qi], 'assay_aligned_lineage_candidate'] = best
    obs.loc[obs.index[qi], 'assay_aligned_lineage_score'] = score
    for lineage, genes in MARKERS.items():
        frac, det, total = fraction_and_detection(data.X, measured, genes)
        obs[f'positive_{lineage.replace(" ", "_")}_fraction'] = frac
        obs[f'positive_{lineage.replace(" ", "_")}_genes'] = det
        obs[f'positive_{lineage.replace(" ", "_")}_counts'] = total
    otherfractions = np.column_stack([obs[f'positive_{lineage.replace(" ", "_")}_fraction'].to_numpy() for lineage in ['Endothelial', 'Fibroblast', 'Pericyte', 'Smooth muscle', 'Immune']])
    ec_regulators = fraction_and_detection(data.X, measured, ['Erg', 'Fli1'])[1]
    ec_vascular = fraction_and_detection(data.X, measured, ['Cdh5', 'Pecam1', 'Kdr', 'Esam', 'Tek', 'Vwf', 'Flt1', 'Tie1'])[1]
    obs['endothelial_regulators_detected'] = ec_regulators
    obs['endothelial_vascular_markers_detected'] = ec_vascular
    # The source paper reports ECs with muscle-associated transcripts. One EC
    # regulator plus three vascular markers warrants review, not a forced CM call.
    obs['CM_EC_program_overlap_review_flag'] = (ec_regulators >= 1) & (ec_vascular >= 3) & obs.positive_CM_fraction.ge(.10).to_numpy()
    cm_positive = obs.positive_CM_genes.ge(3).to_numpy() & obs.positive_CM_counts.ge(10).to_numpy() & obs.positive_CM_fraction.ge(.10).to_numpy() & (obs.positive_CM_fraction.to_numpy() >= 1.5 * otherfractions.max(axis=1)) & keepquery & ~obs.CM_EC_program_overlap_review_flag.to_numpy()
    if reuse_global:
        assert np.array_equal(cm_positive, previous_CM_gate), 'Global reuse requires unchanged measured CM gate.'
    obs['measured_CM_marker_gate'] = cm_positive
    # Positive CM evidence takes precedence over an inconsistent broad vote.
    obs['spatial_lineage_candidate'] = 'Unresolved cell identity'
    obs.loc[~obs.pass_expression_QC, 'spatial_lineage_candidate'] = 'Low-information/QC flag'
    obs.loc[obs.CM_EC_program_overlap_review_flag & obs.pass_expression_QC, 'spatial_lineage_candidate'] = 'CM/EC program overlap; review'
    if reuse_global:
        obs.loc[~cm_positive, 'spatial_lineage_candidate'] = previous_lineage.loc[~cm_positive]
    for lineage in ['Endothelial', 'Fibroblast', 'Pericyte', 'Smooth muscle', 'Immune', 'Epicardial']:
        key = lineage.replace(' ', '_')
        minimum_fraction = .02 if lineage == 'Epicardial' else .05
        consistent = obs[f'positive_{key}_genes'].ge(2).to_numpy() & obs[f'positive_{key}_fraction'].ge(minimum_fraction).to_numpy()
        supported = bool(broadcheck.loc[lineage, 'panel_supported']) if lineage in broadcheck.index else False
        accepted = (best == lineage) & (score >= .65) & (margin >= .15) & consistent[qi] & ~cm_positive[qi] & ~obs.CM_EC_program_overlap_review_flag.to_numpy()[qi] & supported & (not reuse_global)
        obs.loc[obs.index[qi[accepted]], 'spatial_lineage_candidate'] = lineage
    cmref = ref.iloc[ri].broad_class.isin(['Ventricular-enriched CM', 'Atrial-like CM']).to_numpy()
    cmquery = np.flatnonzero(cm_positive[qi])
    obs.loc[obs.index[qi[cmquery]], 'spatial_lineage_candidate'] = 'CM; chamber unresolved'
    if len(cmquery):
        if not reuse_global:
            chamber, chamber_score, chamber_margin, _, _, _ = vote(re[cmref], qe[cmquery], ref.iloc[ri].loc[cmref, 'broad_class'])
            obs['CM_chamber_candidate'] = pd.NA
            obs['CM_chamber_score'] = np.nan
            obs.loc[obs.index[qi[cmquery]], 'CM_chamber_candidate'] = chamber
            obs.loc[obs.index[qi[cmquery]], 'CM_chamber_score'] = chamber_score
        # Chamber candidates need positive measured evidence. Integration alone
        # can call ventricular profiles atrial despite absent atrial transcripts.
        myl2 = obs.measured_Myl2_counts.to_numpy()
        myl4 = obs.measured_Myl4_counts.to_numpy()
        myl7 = obs.measured_Myl7_counts.to_numpy()
        # Myl7 is developmentally expressed in ventricular CMs; never use it
        # as negative evidence against neonatal ventricular identity.
        ventricular_markers = cm_positive & (myl2 >= 2) & (myl2 >= myl4)
        atrial_markers = cm_positive & (myl4 >= 2) & (myl7 >= 1) & (myl2 <= 1)
        obs['measured_ventricular_marker_gate'] = ventricular_markers
        obs['measured_atrial_marker_gate'] = atrial_markers
        obs.loc[ventricular_markers, 'spatial_lineage_candidate'] = 'Ventricular-enriched CM'
        obs.loc[atrial_markers, 'spatial_lineage_candidate'] = 'Atrial-marker-enriched CM'
    np.savez_compressed(CACHE / f'{stage}_broad_alignment_cell_order.npz', reference_cell_ids=ref.iloc[ri].index.to_numpy(dtype=str), query_cell_ids=obs.index[qi].to_numpy(dtype=str))
    print(f'{stage} MARKER_CONFIRMED CM={cm_positive.sum()} chamber assignments {obs.spatial_lineage_candidate.value_counts().to_dict()}', flush=True)
    if only_broad:
        obs.to_csv(folder / 'assay_aligned_broad_diagnostic.tsv.gz', sep='\t', na_rep='NA')
        return
    vref = np.flatnonzero(information & ref.broad_class.eq('Ventricular-enriched CM').to_numpy())
    vquery = np.flatnonzero(obs.spatial_lineage_candidate.eq('Ventricular-enriched CM').to_numpy())
    cr, cq = embedding(corrected[vref], query[vquery], f'{stage}_ventricular_marker_v2_corrected', 20, 30)
    rr, rq = embedding(original[vref], query[vquery], f'{stage}_ventricular_marker_v2_original', 20, 30)
    reference_families = ref.iloc[vref].mapping_family.to_numpy(dtype=str)
    known = reference_families != 'Uncharacterized ventricular candidate'
    famcheck = benchmark(stage, 'family', cr[known], reference_families[known], ref.iloc[vref].loc[known, 'sample'], folder)
    finecheck = benchmark(stage, 'fine_group', cr, ref.iloc[vref].ventricular_group, ref.iloc[vref]['sample'], folder)
    family, fscore, fmargin, _, _, _ = vote(cr[known], cq, reference_families[known])
    rawfamily, _, _, _, _, _ = vote(rr[known], rq, reference_families[known])
    fine, gscore, gmargin, _, _, _ = vote(cr, cq, ref.iloc[vref].ventricular_group)
    rawfine, _, _, _, _, _ = vote(rr, rq, ref.iloc[vref].ventricular_group)
    library_families, library_fine = [], []
    for library in sorted(ref.iloc[vref]['sample'].unique()):
        lib = ref.iloc[vref]['sample'].eq(library).to_numpy()
        value, _, _, _, _, _ = vote(cr[lib & known], cq, reference_families[lib & known])
        library_families.append(value)
        value, _, _, _, _, _ = vote(cr[lib], cq, ref.iloc[vref].loc[lib, 'ventricular_group'])
        library_fine.append(value)
    family_agreement = (library_families[0] == library_families[1]) & (family == library_families[0])
    fine_agreement = (library_fine[0] == library_fine[1]) & (fine == library_fine[0])
    obs['spatial_family_candidate'] = pd.NA
    obs['spatial_family_score'] = np.nan
    obs['spatial_family_raw_reference_agreement'] = False
    obs['spatial_family_library_agreement'] = False
    obs['spatial_mapped_family'] = 'Not assigned ventricular CM'
    obs.loc[obs.index[vquery], 'spatial_mapped_family'] = 'Unresolved ventricular state'
    obs['spatial_fine_group_candidate'] = pd.NA
    obs['spatial_fine_group_score'] = np.nan
    obs['spatial_mapped_fine_group'] = 'Not resolved'
    recover = pd.Series(family).map(famcheck.panel_supported).fillna(False).to_numpy()
    acceptfamily = recover & (fscore >= .60) & (fmargin >= .15) & (family == rawfamily) & family_agreement & obs.iloc[vquery].shared_panel_transcripts.ge(50).to_numpy()
    query_log = normalize(query)
    for program, genes in PROGRAMS.items():
        frac, detected, total = fraction_and_detection(data.X, measured, genes)
        obs[f'measured_{program}_genes'] = detected
        obs[f'measured_{program}_counts'] = total
        obs[f'measured_{program}_fraction'] = frac
        selected = [sharedgenes.index(gene) for gene in genes if gene in sharedgenes]
        score_program = query_log[:, selected].mean(axis=1)
        center = score_program[vquery].mean()
        sd = max(score_program[vquery].std(), 1e-8)
        obs[f'measured_{program}_within_section_vent_z'] = (score_program - center) / sd
    del query_log
    acceptfamily &= (family != 'Cycling-associated') | obs.iloc[vquery].measured_cycling_genes.ge(2).to_numpy()
    acceptfamily &= (family != 'Conduction-program candidate') | obs.iloc[vquery].conduction_panel_genes_detected.ge(2).to_numpy()
    obs.loc[obs.index[vquery], 'spatial_family_candidate'] = family
    obs.loc[obs.index[vquery], 'spatial_family_score'] = fscore
    obs.loc[obs.index[vquery], 'spatial_family_raw_reference_agreement'] = family == rawfamily
    obs.loc[obs.index[vquery], 'spatial_family_library_agreement'] = family_agreement
    obs.loc[obs.index[vquery[acceptfamily]], 'spatial_mapped_family'] = family[acceptfamily]
    catalog = pd.read_csv(OUT / 'reference_group_catalog.tsv', sep='\t').set_index('microenvironment_cluster')
    marker_supported = catalog.loc[fine, 'both_library_raw_corrected_marker_count'].ge(3).to_numpy()
    fine_recover = pd.Series(fine).map(finecheck.panel_supported).fillna(False).to_numpy()
    expected_family = catalog.loc[fine, 'mapping_family'].to_numpy()
    acceptedfine = fine_recover & marker_supported & (gscore >= .70) & (gmargin >= .20) & (fine == rawfine) & fine_agreement & (expected_family == obs.iloc[vquery].spatial_mapped_family.to_numpy()) & obs.iloc[vquery].shared_panel_transcripts.ge(100).to_numpy()
    obs.loc[obs.index[vquery], 'spatial_fine_group_candidate'] = fine
    obs.loc[obs.index[vquery], 'spatial_fine_group_score'] = gscore
    obs.loc[obs.index[vquery[acceptedfine]], 'spatial_mapped_fine_group'] = fine[acceptedfine]
    obs['measured_Ankrd1_program_positive'] = obs.measured_Ankrd1_counts.ge(2) & obs.measured_Ankrd1_associated_genes.ge(3) & obs.measured_Ankrd1_associated_within_section_vent_z.ge(.5) & obs.spatial_lineage_candidate.eq('Ventricular-enriched CM')
    obs['measured_Slit2_program_positive'] = obs.measured_Slit2_counts.ge(1) & obs.measured_Slit2_associated_genes.ge(2) & obs.measured_Slit2_associated_within_section_vent_z.ge(.5) & obs.spatial_lineage_candidate.eq('Ventricular-enriched CM')
    obs['measured_conduction_overlap_positive'] = obs.measured_Cntn2_counts.ge(1) & obs.measured_Gja5_counts.ge(1) & obs.conduction_panel_genes_detected.ge(3) & obs.spatial_lineage_candidate.eq('Ventricular-enriched CM')
    obs['measured_cycling_program_positive'] = obs.measured_cycling_genes.ge(3) & (obs.measured_Mki67_counts.ge(1) | obs.measured_Top2a_counts.ge(1)) & obs.spatial_lineage_candidate.eq('Ventricular-enriched CM')
    obs.to_csv(folder / 'spatial_candidates_and_measured_programs.tsv.gz', sep='\t', na_rep='NA', float_format='%.7g')
    data.obs = obs
    data.obsm['X_assay_aligned_global'] = np.full((len(obs), re.shape[1]), np.nan, dtype=np.float32)
    data.obsm['X_assay_aligned_global'][qi] = qe
    data.obsm['X_assay_aligned_ventricular'] = np.full((len(obs), cr.shape[1]), np.nan, dtype=np.float32)
    data.obsm['X_assay_aligned_ventricular'][vquery] = cq
    data.uns['primary_mapping_columns'] = ['spatial_lineage_candidate', 'spatial_mapped_family', 'spatial_mapped_fine_group']
    data.uns['unintegrated_classifier_columns'] = 'mapped_lineage/mapped_family/mapped_fine_group and associated scores are diagnostic comparisons only; direct transfer showed cross-assay bias.'
    data.write_h5ad(folder / f'{stage}_xenium_counts_and_candidate_mapping.h5ad', compression='gzip')
    summaries = []
    for field in ['spatial_lineage_candidate', 'spatial_mapped_family', 'spatial_mapped_fine_group']:
        for label, n in obs[field].value_counts().items():
            summaries.append({'stage': stage, 'level': field, 'label': label, 'segmented_cells': int(n)})
    pd.DataFrame(summaries).to_csv(folder / 'spatial_candidate_counts.tsv', sep='\t', index=False)
    summary = {'stage': stage, 'source_segmented_cells': len(obs), 'QC_pass': int(obs.pass_expression_QC.sum()), 'measured_CM_marker_gate': int(cm_positive.sum()), 'ventricular_CM_candidates': len(vquery), 'family_assigned': int(acceptfamily.sum()), 'fine_group_assigned': int(acceptedfine.sum()), 'measured_Ankrd1_program_positive': int(obs.measured_Ankrd1_program_positive.sum()), 'measured_Slit2_program_positive': int(obs.measured_Slit2_program_positive.sum()), 'measured_conduction_overlap_positive': int(obs.measured_conduction_overlap_positive.sum()), 'measured_cycling_program_positive': int(obs.measured_cycling_program_positive.sum()), 'runtime_seconds': time.monotonic() - start}
    (folder / 'spatial_summary.json').write_text(json.dumps(summary, indent=2))
    print('SPATIAL_COMPLETE', summary, flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stages', nargs='+', default=['P0', 'P7', 'P14', 'P21'])
    parser.add_argument('--only-broad', action='store_true')
    parser.add_argument('--reuse-global', action='store_true', help='Reuse completed global votes when the measured CM gate and input counts are unchanged.')
    args = parser.parse_args()
    CACHE.mkdir(parents=True, exist_ok=True)
    reference_directory = CACHE if (CACHE / 'reference_shared_panel_counts.npz').exists() else OUT
    cache = np.load(reference_directory / 'reference_shared_panel_counts.npz')
    ref = pd.read_csv(reference_directory / 'reference_metadata.tsv.gz', sep='\t', index_col=0)
    genes = list(cache['genes'])
    for stage in args.stages:
        process(stage, cache, ref, genes, args.only_broad, args.reuse_global)


if __name__ == '__main__':
    main()
