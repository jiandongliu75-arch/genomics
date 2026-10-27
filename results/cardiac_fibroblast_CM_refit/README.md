# Changing fibroblast environments and cardiomyocyte development after birth

This release completes the CM-neighbor refit using the latest fibroblast candidate assignments and directly measured gene patterns. It follows the clear biological style of the preceding report. Original data, spatial coordinates, reviewed CM cohorts, graphs and earlier releases are preserved.

- [Word report](Fibroblast_CM_spatial_reassessment.docx)
- [Full PDF report](Fibroblast_CM_spatial_reassessment.pdf)
- [All six main and one supporting figure plates](Fibroblast_CM_integrated_figures.pdf)
- [Source-table download](Fibroblast_CM_source_tables.zip)
- [Readable report and figure legends](Report_text.md)
- [Numerical validation](numerical_validation.json)
- [Document/release validation](release_validation.json)

## Biological findings

The broader P14 Col8a1 candidates remain more often near arterial structures than matched Gfpt2 candidates, but the contrast decreases to approximately 11% versus 7% from the earlier 19% versus 6%. Most Col8a1 candidates are farther than 50 µm from these structures. This is a vascular-associated component in a broader matrix population.

No fibroblast identity or program association passes the global correction across all 925 primary comparisons. The earlier P7 Postn/cycling-CM metabolic association weakens with broader labels and smaller neighborhoods. P0 elastic-matrix expression and lower working-CM sarcomere expression remain a candidate supported across several sensitivity checks, including CM-gene-excluded scoring; it does not pass the primary correction.

The clearest corrected association is higher sarcomere RNA in P14 Ankrd1 CMs near endocardial-like ECs. It remains with updated FB adjustment and after removing Myl2, but weakens with deposited whole-cell CM RNA. This concerns a component of specialization, not greater overall maturity or proof of contact signaling.

Frozen P7-to-P14 exact neighborhood matches show declining Postn-associated and increasing Col8a1-associated expression across four examined CM/anatomical contexts. Fifteen of 96 sampled neighborhoods contain no eligible FB and the median is two. These descriptions do not track ancestry or establish a fine transition map.

## What was newly fitted

- Candidate, crosschecked and earlier stringent FB identities.
- Seven measured FB patterns: four published marker sets and fibrillar, basement and elastic matrix sets.
- Technical, anatomical and within-region models for separate CM components and the 21-gene index.
- Native/core/depth-matched source checks; native CM outcomes; Myl2 exclusion; CM-gene-excluded rank/control backgrounds; joint-program specificity; same-CM 30/50-µm and shared-boundary comparisons.
- Capillary/endocardial target models with updated FB identity or program adjustment.
- Updated technical-matched FB spatial comparisons and held-out descriptions on frozen exact/Wormhole neighborhood plans.

No new classifier, segmentation, encoder, perfusion/function measurement, biological replicate or lineage experiment was added. Earlier capillary-apposition comparisons are discussed as context and were not rerun here.

## Reading the tables

`fibroblast_assignment_coverage.tsv` uses missing numeric counts and `Not assessable` when a family is absent from an age's reference. These are not biological zeros. Coverage does not establish accuracy; thresholds selected on reference predictions do not yield unbiased prospective spatial performance.

`fibroblast_per_profile_programs.tsv.gz` contains original cell IDs, stage, count view, informative-source flags, depth, detected genes and all seven continuous scores. Scores use matched-control rank excess, not functional activity. Named patterns can overlap. `fibroblast_program_definitions.tsv` gives exact genes. The additional `CM_gene_exclusion_definition.tsv` records the separate source-scoring sensitivity.

`CM_fibroblast_refitted_coefficients.tsv` includes 17,256 estimates, not 17,256 primary discoveries. `primary_analysis=True` identifies the 925-test family; its global BH q values are missing for all other rows. Positive cycling coefficients mean higher cell-cycle RNA. Each exposure is standardized within its age/run; coefficients are not causal effects or a calibrated cross-age developmental scale. Joint-program and geometry sensitivities are separate files. `model_design_diagnostics.tsv` records cohorts, blocks, exposure variation, rank and condition numbers.

Each `{stage}_{geometry}_neighborhood_exposures.npz` stores `cell_ids` and aligned one-dimensional exposures for all frozen CMs, including ineligible ones. Candidate fractions divide by all eligible broad FB neighbors, retaining unassigned cells. Continuous programs average over informative eligible FBs. Coverage/depth/genes are supplied separately. Join by stage and original cell ID, then apply the documented CM cohort and eligibility rules. These arrays are available directly in this folder; the smaller ZIP contains source **tables**, not duplicated NPZ matrices or native counts.

`updated_matched_fibroblast_pairs.tsv.gz` retains paired IDs and reused controls; `updated_spatial_matching_balance.tsv` provides balance. Matched pair count is not independent sample size.

`updated_frozen_transport_neighborhood_descriptors.tsv` updates 96 neighborhoods. `updated_transport_program_comparisons.tsv` reports changes and the supported mass for both saved matching methods. Their underlying RNA information overlaps the molecular representation, so they are not independent validation. Plans and encoder remain in the previous subtype release.

## Reproduce

Use the existing scientific Python with NumPy, pandas, SciPy, anndata, scikit-learn and matplotlib; document Python requires python-docx, Pillow and PyMuPDF, with LibreOffice for rendering. Dependencies also import the preceding analysis definitions; run from the repository root. Local Xenium count inputs and the shared RNA reference must be available at the paths recorded in `input_paths.json`.

```bash
OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 /workspace/heart-atlas/venv/bin/python analysis/run_fibroblast_CM_refit.py --document-python python
```

For an existing completed release, `--resume` verifies source, code and artifact hashes without refitting. The report can be regenerated independently with `python analysis/write_fibroblast_CM_refit_report.py`, LibreOffice conversion, and `python analysis/finalize_fibroblast_CM_refit.py`.

## Interpretation limits and attribution

One available spatial section per age precludes independent-heart inference of developmental effects. Anatomical masks are operational proxies. Cell boundaries and RNA ownership remain uncertain; block intervals do not integrate all segmentation uncertainty. Small marker panels do not establish functional matrix activity, completed CM maturation, causal signaling or cell ancestry.

Original atlas: Wang et al., Nature Cardiovascular Research (2026), https://doi.org/10.1038/s44161-026-00861-z . Perinatal vascular context: Lu et al., Nature Cell Biology (2021), https://doi.org/10.1038/s41556-021-00747-1 . Wormhole: Haviv et al., ICML (2024), https://arxiv.org/abs/2404.09411 . Original source counts and journal PDFs are not redistributed in this new report folder.
