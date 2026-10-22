# Xenium QC: current status and the next segmentation checks

**Completed follow-up:** [Transcript-level boundary inspection and segmentation sensitivity](../xenium_boundary_qc/README.md) now audits the original molecule tables and acquired boundary stains, reconstructs all source matrices, tests boundary alternatives, and supplies targeted review flags while retaining CM profiles without nuclear area. The checks below describe the earlier native-atlas baseline; see the follow-up for the updated status and remaining anatomical limits.

CytoSPACE aligns a single-cell/nucleus reference to spatial expression profiles; it does not repair segmentation. The current whole-heart atlas uses the deposited cell segmentation and count matrices. Its broad identities and spatial associations remain conditional on that segmentation.

## What is implemented

| Check | Current implementation | What it does not establish |
|---|---|---|
| Basic expression information | ≥30 measured transcripts, ≥10 detected genes, positive cell area, finite centroids | A biologically complete or correctly segmented cell |
| Negative-control probes | Control-probe / (gene transcripts + control-probe counts) ≤0.10 | Absence of boundary leakage or an estimate of ambient-RNA contamination |
| Reference information | ≥20 molecules and ≥8 detected genes among the377 fitting targets | Fine-subtype recoverability or correct morphology |
| Reference robustness | Competing broad lineages; weight/margin, original-reference and two-library agreement | Independent cross-platform anatomical accuracy |
| Marker consistency | Lineage anchors; strong CM/non-CM coexpression retained as a review flag | Proof of doublets, or proof that every mixed profile is technical |
| Geometry conservation | Original cell IDs/coordinates/counts conserved; native polygon overlays in selected regions | Correct biological cell borders, complete cells or true nuclei-per-cell counts |
| Available cell-level morphology | Area, nuclear area, nuclear/cell ratio, source control counts, stage/lineage-specific area review | A comprehensive segmentation validation or newly generated boundaries |

[Available-QC report](Xenium_available_QC_diagnostics.pdf), [summary table](available_Xenium_QC_summary.tsv) and [provenance](available_QC_provenance.json) accompany per-stage `available_Xenium_QC_metrics.tsv.gz` files. A [stratified candidate review sample](candidate_segmentation_review_sample.tsv.gz) contains original cell IDs and locations for a future visual audit; these profiles have **not** been manually adjudicated.

The prior basic filters are operational analysis thresholds, not universal Xenium acceptance standards. All original cells remain available, including failed/uncertain cells. Missing nuclear area is flagged rather than automatically excluded. It occurs in approximately 10.5%, 14.9%, 11.4% and 15.7% of the P0/P7/P14/P21 profiles. Section truncation, cytoplasmic profiles and segmentation can all contribute; these observations do not diagnose the cause.

The measured-marker plots also show CM-marker RNA in multiple non-CM classes after conservative mixed-profile flags. Detection alone does not establish contamination, and high-abundance CM genes need quantitative interpretation. This residual signal makes cell-interface transcript localization a particularly relevant next check before fine non-CM states or physical cell-contact claims are accepted. The audit does not tune thresholds to force visually pure marker plots.

## Recommended next checks, in priority order

1. **Audit segmentation where conclusions depend on cell boundaries.** Overlay nuclear and cell polygons, the highest-resolution available morphology/membrane images, and individual transcript coordinates. Sample ordinary and flagged regions across ages and lineages, including CM–EC, CM–fibroblast and peri-arterial interfaces. Score merged neighboring cells, over-split CMs, truncated profiles, misplaced boundaries, excluded nuclei and cross-boundary transcript assignment. A DAPI-only overlay is useful for nuclei but cannot validate the full cytoplasmic outline of a cardiomyocyte.
2. **Review decoding, controls and spatial artifacts.** Examine run-level decoding and negative-probe/negative-codeword metrics, transcript QV distributions, assay coverage, imaging focus, tissue folds/tears, autofluorescence, acquisition-tile seams and transcript/morphology registration. Confirm the quality cutoff used to generate the deposited matrix from run metadata and molecule-level tables; Xenium commonly uses QV≥20 for high-quality transcripts, but that has not been independently verified here. Check regional patterns rather than relying only on global averages.
3. **Quantify molecule assignment at cell interfaces.** From molecule coordinates and assignments, measure unassigned transcript fractions, transcript-to-nucleus and transcript-to-boundary distances, and lineage-marker distributions in inner versus outer cell regions. Look for small EC or fibroblast profiles receiving disproportionately peripheral CM transcripts. Nuclear/cytoplasmic localization and genuine shared expression require gene- and cell-type-aware interpretation. The `unassigned_codeword_counts` column in cell summaries is **not** the fraction of transcripts unassigned to cells.
4. **Use cell-type-aware cell QC.** Inspect transcript counts, detected genes, area, nuclear/cell area ratio and transcript density by stage and broad lineage. Use robust review flags, not an indiscriminate global upper-count or upper-size cutoff. The targeted gene panel does not justify importing standard scRNA mitochondrial-percentage cutoffs without checking which mitochondrial genes are actually measured.
5. **Compare segmentation/assignment alternatives in representative regions.** First assess improved vendor segmentation if compatible morphology channels and raw inputs exist. An independent morphology- or transcript-aware method can then provide a sensitivity comparison, with parameters suited to cardiac tissue. The goal is not simply higher marker purity: assess morphology, counts retained, merged/split objects and robust biological assignments together. Choose settings without maximizing the desired Ankrd1 or subtype association.
6. **Repeat consequential biology after the audit.** Rebuild counts from accepted alternative segmentations/assignments, repeat broad identities, and retest subtype placements and the Ankrd1 peri-arterial comparison. Also compare conservative subsets, such as non-mixed profiles with strong marker support. For CMs, a nucleus-centered subset is a sensitivity analysis, not necessarily the most complete representation of the cell.

## Cardiomyocyte-specific cautions

- Postnatal CMs are large and increasingly binucleated/multinucleated. Multiple nuclei, high RNA counts or large area alone do not establish a doublet.
- A CM can be split into several spatial segments, or a nearby EC can receive CM cytoplasmic RNA. Segmented objects are not automatically a one-to-one biological cell census.
- Strong CM/non-CM coexpression is a review trigger. Droplet-based Scrublet/scDblFinder calls are not a substitute for examining spatial geometry.
- snRNA ambient-RNA correction should not be applied indiscriminately to Xenium counts. Investigate decoding background and spatial transcript misassignment first. If reassignment/correction is justified, preserve original counts and report sensitivity to the correction.

Molecule-level QV, transcript assignment, segmentation alternatives and an independently scored morphology audit are still outstanding. The present maps and available-QC diagnostics are reviewable starting points, not evidence that these remaining checks have already been completed.
