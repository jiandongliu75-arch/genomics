"""Same-stage, shared-panel reference transfer with explicit unresolved labels.

Expression-only mapping; no spatial smoothing, no Kcnd2 imputation, and no
assumption that reference softmax scores are calibrated Xenium probabilities.
"""
from pathlib import Path
import argparse
import json
import time
import warnings

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import classification_report, confusion_matrix
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.exceptions import ConvergenceWarning

ROOT = Path(__file__).resolve().parents[1]
BASE = Path('/workspace/heart-atlas')
OUT = ROOT / 'results/xenium_ventricular_mapping'
CACHE = BASE / 'results/xenium_mapping_cache'
CM_GENES = 'Tnnt2 Tnni3 Myh6 Myh7 Myl2 Actc1'.split()
CYCLING_GENES = 'Mki67 Top2a Cenpf Cenpe Cdca2 Kif11 Kif15 Smc4 Hjurp Atad2'.split()
CONDUCTION_GENES = 'Slit2 Sema3a Ephb1 Gja5 Cntn2 Robo1'.split()


def normalize(counts):
    counts = np.asarray(counts, dtype=np.float32)
    total = counts.sum(axis=1)
    return np.log1p(counts * (10000 / np.maximum(total, 1))[:, None])


def fit_model(expression, labels):
    model = Pipeline([('scale', StandardScaler()), ('classifier', LogisticRegression(C=.2, class_weight='balanced', solver='lbfgs', max_iter=1500, tol=1e-4))])
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always', ConvergenceWarning)
        model.fit(expression, labels)
    if any(issubclass(item.category, ConvergenceWarning) for item in caught):
        model.set_params(classifier__max_iter=4000, classifier__warm_start=True)
        model.fit(expression, labels)
        if model.named_steps['classifier'].n_iter_.max() >= 4000:
            raise RuntimeError('Classifier did not converge; do not publish this fit.')
    return model


def thin(counts, seed):
    rng = np.random.default_rng(seed)
    totals = counts.sum(axis=1)
    result = np.zeros(counts.shape, dtype=np.float32)
    for i in range(len(counts)):
        if totals[i] > 0:
            n = min(100, max(1, int(round(float(totals[i])))))
            probabilities = counts[i].astype(float) / totals[i]
            probabilities /= probabilities.sum()
            result[i] = rng.multinomial(n, probabilities)
    return normalize(result)


def marker_count(matrix, genes, wanted):
    ids = [genes.index(gene) for gene in wanted if gene in genes]
    return np.asarray((matrix[:, ids] >= 1).sum(axis=1)).ravel()


def evaluate_and_fit(stage, level, counts, original, labels, libraries, folder):
    labels = np.asarray(labels, dtype=str)
    libraries = np.asarray(libraries, dtype=str)
    classes = np.unique(labels)
    expression = normalize(counts)
    low = thin(counts, 42 + int(stage[1:]))
    rows = []
    confusions = []
    for heldout in sorted(np.unique(libraries)):
        train = libraries != heldout
        test = ~train
        print(f'{stage} {level} HELDOUT {heldout} train={train.sum()} test={test.sum()}', flush=True)
        model = fit_model(expression[train], labels[train])
        for condition, values in [('full_shared_panel', expression[test]), ('multinomial_up_to_100_molecules', low[test])]:
            pred = model.predict(values)
            report = classification_report(labels[test], pred, labels=classes, output_dict=True, zero_division=0)
            for label in classes:
                item = report[label]
                rows.append({'stage': stage, 'level': level, 'heldout_library': heldout, 'condition': condition, 'reference_class': label, 'precision': item['precision'], 'recall': item['recall'], 'f1': item['f1-score'], 'test_nuclei': int(item['support']), 'training_nuclei': int((labels[train] == label).sum())})
            cm = confusion_matrix(labels[test], pred, labels=classes)
            for i, truth in enumerate(classes):
                for j, prediction in enumerate(classes):
                    confusions.append({'stage': stage, 'level': level, 'heldout_library': heldout, 'condition': condition, 'truth': truth, 'predicted': prediction, 'nuclei': int(cm[i, j])})
    metrics = pd.DataFrame(rows)
    metrics.to_csv(folder / f'{level}_heldout_metrics.tsv', sep='\t', index=False)
    pd.DataFrame(confusions).to_csv(folder / f'{level}_heldout_confusion.tsv', sep='\t', index=False)
    summary = metrics.groupby('reference_class').agg(min_precision=('precision', 'min'), min_recall=('recall', 'min'), min_f1=('f1', 'min'), min_test_nuclei=('test_nuclei', 'min'), min_training_nuclei=('training_nuclei', 'min'))
    summary['panel_recoverable'] = summary.min_precision.ge(.60) & summary.min_recall.ge(.30) & summary.min_f1.ge(.40) & summary.min_test_nuclei.ge(30) & summary.min_training_nuclei.ge(30)
    summary.to_csv(folder / f'{level}_panel_recoverability.tsv', sep='\t')
    print(f'{stage} {level} RECOVERABLE {summary.index[summary.panel_recoverable].tolist()}', flush=True)
    model = fit_model(expression, labels)
    rawmodel = fit_model(normalize(original), labels)
    for kind, fitted in [('corrected_reference', model), ('original_reference', rawmodel)]:
        classifier, scale = fitted.named_steps['classifier'], fitted.named_steps['scale']
        coefficient = classifier.coef_ / scale.scale_[None, :]
        intercept = classifier.intercept_ - coefficient @ scale.mean_
        np.savez_compressed(folder / f'{level}_{kind}_coefficients.npz', classes=fitted.classes_, coefficient=coefficient, intercept=intercept, fit_iterations=classifier.n_iter_)
    return model, rawmodel, summary


def prediction(model, rawmodel, expression):
    scores = model.predict_proba(expression)
    order = np.argsort(scores, axis=1)
    first, second = order[:, -1], order[:, -2]
    names = model.classes_[first]
    top = scores[np.arange(len(scores)), first]
    margin = top - scores[np.arange(len(scores)), second]
    raw = rawmodel.predict(expression)
    return names, top, margin, raw, scores


def process(stage, cache, reference, genes, catalog):
    start = time.monotonic()
    folder = OUT / stage
    folder.mkdir(exist_ok=True)
    mask = reference.age.eq(stage).to_numpy()
    corrected, original = cache['corrected'][mask], cache['original'][mask]
    ref = reference.loc[mask].copy()
    info = (corrected.sum(axis=1) >= 10) & ((corrected >= 1).sum(axis=1) >= 5)
    broadmask = info & ref.broad_class.ne('Unresolved/mixed').to_numpy()
    broadmodel, broadraw, broadcheck = evaluate_and_fit(stage, 'lineage', corrected[broadmask], original[broadmask], ref.loc[broadmask, 'broad_class'], ref.loc[broadmask, 'sample'], folder)
    ventmask = info & ref.broad_class.eq('Ventricular-enriched CM').to_numpy()
    fammask = ventmask & ref.mapping_family.ne('Uncharacterized ventricular candidate').to_numpy()
    fammodel, famraw, famcheck = evaluate_and_fit(stage, 'family', corrected[fammask], original[fammask], ref.loc[fammask, 'mapping_family'], ref.loc[fammask, 'sample'], folder)
    finemodel, fineraw, finecheck = evaluate_and_fit(stage, 'fine_group', corrected[ventmask], original[ventmask], ref.loc[ventmask, 'ventricular_group'], ref.loc[ventmask, 'sample'], folder)
    groupcat = catalog.set_index('microenvironment_cluster')
    finecheck['upstream_both_library_marker_count'] = groupcat.loc[finecheck.index, 'both_library_raw_corrected_marker_count']
    finecheck['upstream_support_status'] = groupcat.loc[finecheck.index, 'support_status']
    finecheck['panel_recoverable'] &= finecheck.upstream_both_library_marker_count.ge(3)
    finecheck.to_csv(folder / 'fine_group_panel_recoverability.tsv', sep='\t')
    stage_pattern = 'WT_P07__*' if stage == 'P7' else f'WT_{stage}__*'
    source = next((BASE / 'data/xenium').glob(stage_pattern))
    with h5py.File(source / 'cell_feature_matrix.h5') as handle:
        m = handle['matrix']
        matrix = sparse.csc_matrix((m['data'][:], m['indices'][:], m['indptr'][:]), shape=tuple(m['shape'][:])).T.tocsr()
        barcodes = np.array([v.decode() for v in m['barcodes'][:]])
        names = np.array([v.decode() for v in m['features/name'][:]])
        feature_ids = np.array([v.decode() for v in m['features/id'][:]])
        types = np.array([v.decode() for v in m['features/feature_type'][:]])
    cells = pd.read_csv(source / 'cells.csv.gz').set_index('cell_id').loc[barcodes].copy()
    gene_mask = types == 'Gene Expression'
    gene_counts = matrix[:, gene_mask].tocsr()
    measured = list(names[gene_mask])
    shared = matrix[:, [list(names).index(gene) for gene in genes]].toarray().astype(np.float32)
    expression = normalize(shared)
    assert np.array_equal(np.asarray(gene_counts.sum(axis=1)).ravel(), cells.transcript_counts.to_numpy())
    cells['stage'] = stage
    cells['section'] = source.name
    cells['measured_gene_count'] = gene_counts.getnnz(axis=1)
    cells['shared_panel_transcripts'] = shared.sum(axis=1)
    cells['shared_panel_detected_genes'] = (shared >= 1).sum(axis=1)
    cells['negative_probe_fraction'] = cells.control_probe_counts / np.maximum(cells.transcript_counts + cells.control_probe_counts, 1)
    cells['missing_nucleus_area_flag'] = cells.nucleus_area.isna()
    cells['valid_coordinates'] = np.isfinite(cells.x_centroid) & np.isfinite(cells.y_centroid)
    cells['pass_expression_QC'] = cells.transcript_counts.ge(30) & cells.measured_gene_count.ge(10) & cells.cell_area.gt(0) & cells.valid_coordinates & cells.negative_probe_fraction.le(.1)
    cells['pass_transfer_information'] = cells.pass_expression_QC & cells.shared_panel_transcripts.ge(20) & cells.shared_panel_detected_genes.ge(8)
    cells['CM_markers_detected'] = marker_count(gene_counts, measured, CM_GENES)
    cells['cycling_markers_detected'] = marker_count(gene_counts, measured, CYCLING_GENES)
    cells['conduction_panel_genes_detected'] = marker_count(gene_counts, measured, CONDUCTION_GENES)
    fibmarkers = marker_count(gene_counts, measured, 'Col1a1 Col1a2 Col3a1 Dcn Lum Pdgfra Tcf21'.split())
    endomarkers = marker_count(gene_counts, measured, 'Pecam1 Cdh5 Kdr Esam Tek'.split())
    def marker_fraction(wanted):
        ids = [measured.index(gene) for gene in wanted if gene in measured]
        return np.asarray(gene_counts[:, ids].sum(axis=1)).ravel() / np.maximum(cells.transcript_counts.to_numpy(), 1)
    cells['CM_marker_transcript_fraction'] = marker_fraction(CM_GENES)
    cells['fibroblast_marker_transcript_fraction'] = marker_fraction('Col1a1 Col1a2 Col3a1 Dcn Lum Pdgfra Tcf21'.split())
    cells['endothelial_marker_transcript_fraction'] = marker_fraction('Pecam1 Cdh5 Kdr Esam Tek'.split())
    cells['CM_nonCM_marker_codetection_flag'] = cells.CM_markers_detected.ge(2) & ((fibmarkers >= 3) | (endomarkers >= 3))
    cells['CM_and_nonCM_marker_overlap_flag'] = cells.CM_markers_detected.ge(2) & cells.CM_marker_transcript_fraction.ge(.15) & (((fibmarkers >= 3) & cells.fibroblast_marker_transcript_fraction.ge(.15)) | ((endomarkers >= 3) & cells.endothelial_marker_transcript_fraction.ge(.15)))
    best, score, margin, rawbest, _ = prediction(broadmodel, broadraw, expression)
    cells['lineage_candidate'] = best
    cells['lineage_score'] = score
    cells['lineage_margin'] = margin
    cells['lineage_raw_reference_agreement'] = best == rawbest
    recover = pd.Series(best).map(broadcheck.panel_recoverable).fillna(False).to_numpy()
    accepted = cells.pass_transfer_information.to_numpy() & recover & (score >= .70) & (margin >= .20) & (best == rawbest)
    iscm = np.isin(best, ['Ventricular-enriched CM', 'Atrial-like CM'])
    accepted &= ~iscm | cells.CM_markers_detected.ge(2).to_numpy()
    cells['mapped_lineage'] = np.where(accepted, best, 'Unresolved cell identity')
    cells.loc[~cells.pass_expression_QC, 'mapped_lineage'] = 'Low-information/QC flag'
    ventquery = cells.mapped_lineage.eq('Ventricular-enriched CM').to_numpy()
    cells['family_candidate'] = pd.NA
    cells['family_score'] = np.nan
    cells['family_margin'] = np.nan
    cells['family_raw_reference_agreement'] = False
    cells['family_panel_recoverable'] = False
    cells['mapped_family'] = 'Not assigned ventricular CM'
    cells.loc[ventquery, 'mapped_family'] = 'Unresolved ventricular state'
    cells['fine_group_candidate'] = pd.NA
    cells['fine_group_score'] = np.nan
    cells['fine_group_margin'] = np.nan
    cells['fine_group_raw_reference_agreement'] = False
    cells['fine_group_panel_recoverable'] = False
    cells['mapped_fine_group'] = 'Not resolved'
    if ventquery.any():
        idx = cells.index[ventquery]
        best, score, margin, rawbest, _ = prediction(fammodel, famraw, expression[ventquery])
        recover = pd.Series(best).map(famcheck.panel_recoverable).fillna(False).to_numpy()
        cells.loc[idx, 'family_candidate'] = best
        cells.loc[idx, 'family_score'] = score
        cells.loc[idx, 'family_margin'] = margin
        cells.loc[idx, 'family_raw_reference_agreement'] = best == rawbest
        cells.loc[idx, 'family_panel_recoverable'] = recover
        acceptfamily = recover & (score >= .60) & (margin >= .15) & (best == rawbest)
        # Conservative measured-marker consistency gates; no imputed evidence.
        acceptfamily &= (best != 'Cycling-associated') | cells.loc[idx, 'cycling_markers_detected'].ge(2).to_numpy()
        acceptfamily &= (best != 'Conduction-program candidate') | cells.loc[idx, 'conduction_panel_genes_detected'].ge(2).to_numpy()
        cells.loc[idx[acceptfamily], 'mapped_family'] = best[acceptfamily]
        best, score, margin, rawbest, _ = prediction(finemodel, fineraw, expression[ventquery])
        recover = pd.Series(best).map(finecheck.panel_recoverable).fillna(False).to_numpy()
        cells.loc[idx, 'fine_group_candidate'] = best
        cells.loc[idx, 'fine_group_score'] = score
        cells.loc[idx, 'fine_group_margin'] = margin
        cells.loc[idx, 'fine_group_raw_reference_agreement'] = best == rawbest
        cells.loc[idx, 'fine_group_panel_recoverable'] = recover
        group_family = groupcat.loc[best, 'mapping_family'].to_numpy()
        acceptfine = recover & (score >= .65) & (margin >= .15) & (best == rawbest) & (group_family == cells.loc[idx, 'mapped_family'].to_numpy())
        cells.loc[idx[acceptfine], 'mapped_fine_group'] = best[acceptfine]
    for gene in 'Ankrd1 Slit2 Gja5 Cntn2 Sema3a Ephb1 Ddc Ncam1 Pola1 Rad51b Mki67 Top2a Tnnt2 Tnni3 Myl2 Myl4 Myl7 Npr3 Upk3b Wt1 Thsd4 Ppara Rreb1'.split():
        if gene in measured:
            cells[f'measured_{gene}_counts'] = gene_counts[:, measured.index(gene)].toarray().ravel()
    cells.to_csv(folder / 'cell_mapping_and_QC.tsv.gz', sep='\t', na_rep='NA', float_format='%.7g')
    var = pd.DataFrame({'gene_symbol': np.array(measured), 'source_feature_id': feature_ids[gene_mask], 'exact_shared_transfer_feature': [gene in genes for gene in measured]}, index=np.array(measured))
    mapped = ad.AnnData(gene_counts, obs=cells, var=var)
    mapped.obsm['spatial'] = cells[['x_centroid', 'y_centroid']].to_numpy()
    mapped.uns['spatial_coordinate_units'] = 'micrometers; native Xenium coordinates; y increases down'
    mapped.uns['Kcnd2_status'] = 'NOT ASSAYED. No Kcnd2 measurement or imputation present.'
    mapped.uns['mapping_limits'] = 'Reference-derived candidate labels. Scores not calibrated Xenium probabilities. Fine labels withheld unless panel/score/raw-reference/family criteria pass. No spatial smoothing.'
    mapped.write_h5ad(folder / f'{stage}_xenium_counts_and_candidate_mapping.h5ad', compression='gzip')
    summaries = []
    for field in ['mapped_lineage', 'mapped_family', 'mapped_fine_group']:
        for label, n in cells[field].value_counts().items():
            summaries.append({'stage': stage, 'level': field, 'label': label, 'segmented_cells': int(n)})
    pd.DataFrame(summaries).to_csv(folder / 'mapping_counts.tsv', sep='\t', index=False)
    summary = {'stage': stage, 'source_segmented_cells': len(cells), 'pass_expression_QC': int(cells.pass_expression_QC.sum()), 'pass_transfer_information': int(cells.pass_transfer_information.sum()), 'ventricular_CM_candidates': int(ventquery.sum()), 'ventricular_family_assigned': int((ventquery & ~cells.mapped_family.eq('Unresolved ventricular state').to_numpy()).sum()), 'fine_group_assigned': int(cells.mapped_fine_group.ne('Not resolved').sum()), 'missing_nucleus_area_not_automatically_excluded': int(cells.missing_nucleus_area_flag.sum()), 'median_transcripts': float(cells.transcript_counts.median()), 'median_measured_genes': float(cells.measured_gene_count.median()), 'CM_nonCM_marker_overlap_flag': int(cells.CM_and_nonCM_marker_overlap_flag.sum()), 'runtime_seconds': time.monotonic() - start}
    (folder / 'section_summary.json').write_text(json.dumps(summary, indent=2))
    print('SECTION_COMPLETE', summary, flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stages', nargs='+', default=['P0', 'P7', 'P14', 'P21'])
    args = parser.parse_args()
    cache = np.load(CACHE / 'reference_shared_panel_counts.npz')
    ref = pd.read_csv(CACHE / 'reference_metadata.tsv.gz', sep='\t', index_col=0)
    assert np.array_equal(cache['cell_ids'], ref.index)
    genes = list(cache['genes'])
    catalog = pd.read_csv(OUT / 'reference_group_catalog.tsv', sep='\t')
    for stage in args.stages:
        process(stage, cache, ref, genes, catalog)
    settings = {'classifier': 'Same-stage multinomial logistic regression; exact shared targets only; reference-fitted gene standardization, class_weight=balanced, C=0.2, lbfgs, max_iter=1500 (retry4000); convergence required. Coefficient exports absorb reference means/scales and act directly on log-normalized shared-panel inputs.', 'benchmark': 'Hold out each of the two original same-age libraries in turn. Test full corrected shared-panel profiles and a multinomial technical simulation with up to 100 molecules per reference nucleus; no upsampling of low-count nuclei. Libraries not assumed independent animals. This is reference-panel recoverability, not cross-platform spatial accuracy.', 'recoverability': 'Minimum precision>=0.60, recall>=0.30, F1>=0.40 across both libraries and both coverage conditions; >=30 train/test nuclei per class. Fine groups additionally require >=3 upstream both-library raw/corrected markers.', 'query_QC': 'Measured transcripts>=30, measured genes>=10, positive cell area, finite centroids, negative-probe fraction<=0.10. Transfer additionally shared transcripts>=20 and shared genes>=8. Missing nuclear area is flagged, not automatically excluded in a thin section.', 'lineage_acceptance': 'Panel recoverability, score>=0.70, margin>=0.20, agreement between corrected-reference and original-reference fits. CM assignments also require >=2 distinct measured CM markers.', 'family_acceptance': 'Only accepted ventricular CMs; recoverability, score>=0.60, margin>=0.15, original-reference agreement; cycling/conduction candidates also >=2 corresponding measured genes.', 'fine_acceptance': 'Recoverability plus upstream marker count, score>=0.65, margin>=0.15, original-reference agreement and consistency with accepted family.', 'interpretation': 'Scores are uncalibrated normalized classifier outputs, not Xenium posterior probabilities. Labels are candidates. Family merging reduces unsupported fine distinctions.', 'spatial': 'No spatial smoothing or neighbor-voting labels. Native centroid coordinates; acquired morphology provides image context only.', 'assay_limits': 'Kcnd2, Irx5 and other absent genes are not measured, imputed, or used to orient anatomy. Original integer Xenium counts unchanged. Controls/QC handled as imaging assay; droplet doublet/ambient methods not applied blindly.', 'segmentation_limits': 'Simple CM/non-CM co-detection is descriptive. Strong overlap flag additionally requires >=15% of measured transcripts in CM markers and >=15% in a fibroblast/endothelial marker set with >=3 corresponding genes. Neither establishes doublets; real shared programs and boundary leakage are alternatives. No transcript-level reassignment or new segmentation performed.', 'stage_limits': 'One deposited acquisition region per stage here; section/animal replication not established. Spatial counts are not an unbiased abundance census or independent animals.'}
    (OUT / 'mapping_methods.json').write_text(json.dumps(settings, indent=2))


if __name__ == '__main__':
    main()
