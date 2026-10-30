# Alternative fibroblast identities and overlapping programs

This focused reassessment compares three different RNA-reference classifiers with direct scoring of the atlas paper's measured marker modules. It assigns candidate transcriptional states and locations, not lineage-traced developmental fate. Earlier atlas labels, segmentation results, CM maturation coefficients and transport plans remain unchanged.

## Read and download

- [Editable Word report](Fibroblast_alternative_identity_reassessment.docx)
- [Rendered report PDF](Fibroblast_alternative_identity_reassessment.pdf)
- [Four integrated figure plates](Fibroblast_multimethod_integrated_figures.pdf)
- [Report text](Report_text.md)
- [Tables, per-cell annotations, programs and model arrays](Fibroblast_multimethod_source_data.zip)
- [Publication figures: PNG, vector PDF and editable SVG](Fibroblast_multimethod_publication_figures.zip)

## What changes biologically

At P21, the measured Gfpt2/Pi16/Adgrd1 marker program is enriched in **1,339 of 10,867 eligible FB profiles**, including **327** also enriched for the Col8a1 marker program. Among 2,716 eligible public Gfpt2-labelled profiles, **1,058** meet the Gfpt2-program flag. These are descriptive three-/four-gene enrichment rules, not validated cell-type calls or estimates of biological abundance.

A narrower Gfpt2 candidate subset is supported by reference models: **85** eligible P21 profiles pass the multi-method cell screen without a supported core conflict, and **75** retain cell/core support and the same majority label after naming genes are withheld. These numbers measure assignment support, not the total true Gfpt2 population.

The 280 P21 profiles called Gfpt2 by the naming-marker-withheld majority still show approximately 84% Gfpt2 and 91% Pi16 detection when inspected afterward. These genes did not enter those grouping models. The RNA reference also contains Gfpt2 expression within broader Col8a1 and Lsamp groups, and the publication lists Col8a1 among its Lsamp markers. A narrow family and a broader overlapping program need not coincide.

New models recover broader Col8a1- and Lsamp-associated candidate maps. They do not prove every new assignment correct. P21 Gfpt2 precision remains limited even in opposite-library RNA tests; selecting whichever model calls the most Gfpt2 cells would overstate resolution.

## Approaches

| Approach | Evidence | Main qualification |
|---|---|---|
| Balanced multinomial logistic classifier | Supervised multi-gene separation with equal class weighting | Training-reference transform can retain platform mismatch |
| Rank-centroid classifier | Relative expression within each cell and family rank templates | High recall can accompany false-positive Gfpt2 calls |
| Negative-binomial count classifier | Library-balanced expression templates and measured counts | Count-distribution/reference assumptions remain |
| Published marker-module scoring | Direct measured signatures independent of our family labels | Small, overlapping modules measure program evidence rather than fate |

Unlike the earlier projection, the new classifiers do not separately recenter query genes or require a naming marker in each cell. Models use training-reference transformations and do not impose a target spatial family proportion. Same-age references are primary; P14/P21 pooled references, naming-marker withholding, raw RNA references, empirical-RNA-prior proposals and information-loss tests are sensitivities.

## Candidate criteria and uncertainty

- Source eligibility retains frozen broad-FB segmentation/marker QC. All broad FB profiles remain in annotation files.
- A cell candidate requires two of three models to pass operational family-specific reference confidence screens and agree, with sufficient measured gene information. The screen is selected from held-out RNA predictions; it is not calibrated Xenium accuracy.
- A supported core majority for another family triggers review. Cell/core and naming-marker-withheld agreement define stricter subsets.
- Nuclear support is optional. Supported nuclear alternatives remain explicit review fields; they are not silently erased. Even a crosschecked candidate can retain this review flag, which should be tested as an exclusion in downstream sensitivities.
- Public labels and marker-module scores remain separately accessible. No forced family quota, spatial smoothing or anatomical valve label is applied.

Use family candidates for questions about discrete RNA-state identity and continuous module scores for programs shared across states. The prior CM association models have not been recomputed using these new layers.

## Data layout

`{stage}_multimethod_FB_annotations.tsv.gz` contains original profile IDs/coordinates, all method/view/variant proposals, candidate layers and review flags. `{stage}_{native|refined|core|nuclear}_published_module_scores.tsv.gz` contains direct program scores, control percentiles and measured-marker counts. `{stage}_published_module_scores.tsv.gz` is the identical refined-view convenience export.

`{stage}_model_predictions.npz` stores labelled axes: **variant × view × method × cell**. `families[label_code]` recovers categorical proposals. Confidence/margin arrays are float32; support/eligibility flags and categorical calls are retained. Example:

```python
import numpy as np
d = np.load('P21_model_predictions.npz')
labels = d['families'][d['label_code'][0, 1, 0]]
cell_ids = d['cell_ids']
print(d['variants'], d['views'], d['methods'])
```

Tables include [reference benchmarks](RNA_opposite_library_model_benchmark.tsv), [training genes](RNA_fold_training_genes.tsv), [candidate coverage](final_candidate_coverage.tsv), [withheld-gene expression](heldout_naming_marker_profiles.tsv), [RNA subcluster marker overlap](RNA_reference_subcluster_marker_audit.tsv), [program overlap](published_Gfpt2_Col8a1_program_overlap.tsv), [public label comparisons](public_Gfpt2_assignment_comparison.tsv) and [boundary versus depth-loss stability](core_vs_depth_matched_assignment_stability.tsv).

The publication and reviewed RNA counts involve different QC/annotation cohorts; the [count comparison](reviewed_vs_published_RNA_family_counts.tsv) is not a matched-cell reclassification table. A source-author RNA barcode annotation crosswalk was not recovered.

## Reproduce

```bash
OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 \
  /workspace/heart-atlas/venv/bin/python analysis/run_fibroblast_multimethod_identity.py --document-python python
```

The scientific environment needs NumPy, pandas, SciPy, scikit-learn, AnnData and matplotlib. Document Python needs python-docx, Pillow and PyMuPDF; LibreOffice renders layout checks. Validation independently reconstructs categorical proposals/consensus, recalculates reference metrics, checks omitted genes and coordinates, compares NB inference with an independent log-PMF calculation, and tests query-composition invariance. Reference recovery and computational agreement do not create independent Xenium truth or biological replication. One available section per age remains a material limitation.

## Attribution

Original atlas: Wang et al., *Spatially guided in vivo single-cell functional genomics of postnatal heart*, Nature Cardiovascular Research (2026), https://doi.org/10.1038/s44161-026-00861-z . Published marker/cell-count definitions are from Supplementary Tables 1–3. Public data derivatives retain attribution. Journal PDFs are not redistributed.
