# Reconstructed atlas relationships and anatomy-led cardiac development analysis

[Eight-page report](cardiac_spatial_reconstruction_report.pdf) · [Selected experiment](EXPERIMENT_PLAN.md) · [Exact-replication limits](AUTHOR_INPUTS_NOT_PUBLIC.md)

The computational roadmap is implemented: reconstruct the released definitions, recover the descriptive relationship, test shared anatomy, measure vascular access and apposition separately, examine different CM programs, infer conditional intermediate RNA programs, and select an experiment. **Exact original Figure 4e replication remains impossible from the released inputs:** the anatomical `source` assignments, precomputed maturation score and neighbor-level table are absent. No inferred day is treated as an observed sample or a new biological replicate.

## Findings that change the biological question

**The atlas's descriptive direction is recovered.** In the independent ventricular selector at P7, more distant capillary EC, proliferating EC and pericyte neighborhoods have lower measured CM maturation-index values. For capillary ECs in the explicit 6NN reconstruction, cell-level Spearman is approximately −0.074, whereas the Pearson correlation of distance-bin means is approximately −0.883. A visually strong binned curve therefore coexists with modest individual-cell association. The tissue-screened Delaunay alternative gives the same direction.

**A large part of that pattern is regional.** Approximately 71% of the P7 capillary-area/maturation-index covariance lies between 200µm regions. This is an exact descriptive covariance partition, not 71% of a causal effect or variance explained. We then compare cells within local regions and use explicit technical, morphology and CM-state sensitivities.

**More apposition does not associate uniformly with every maturation component.** After matching the same local 200µm region and author CM state, restricting separation to ≤100µm, and tightly matching capture, CM area/perimeter and nearby capillary area, the highlighted core-RNA contrasts are:

| Stage | Matched pairs | Metabolic-regulatory RNA difference, SD [95% block interval] | Lower-cycling RNA difference, SD [95% block interval] |
|---|---:|---:|---:|
| P0 | 849 | +0.000 [−0.072, +0.071] | +0.011 [−0.053, +0.082] |
| P7 | 1,291 | +0.037 [−0.013, +0.101] | +0.062 [+0.011, +0.121] |
| P14 | 1,459 | +0.109 [+0.050, +0.166] | +0.033 [−0.019, +0.084] |
| P21 | 2,313 | +0.086 [+0.033, +0.143] | −0.018 [−0.071, +0.035] |

These contrasts compare high with zero sampled capillary apposition. They are exploratory within-section associations; the displayed intervals are not multiplicity-adjusted animal-level confirmation. Full structural, excitation/contraction, native/core results are retained. After the strict match, aggregate standardized imbalances are below 0.016 for all four matching variables; nearby capillary-area imbalance is below 0.004. The initial looser match showed residual capillary-area imbalance and was replaced before release.

**Geometry and identities matter.** Conditional metabolic associations at later ages are positive under native polygon and refined-boundary definitions, but effect magnitudes differ. Strict nuclear-supported capillary definitions attenuate the associations, especially at P7/P14. These identity filters can remove both ambiguous profiles and potentially integrated EC states; no definition is treated as ground truth.

**Nearby EC and CM cycling is another, distinct early relationship.** Nuclear-clean EC sources and nuclear CM RNA show small positive cycling associations within 50µm at P0/P7 (approximately +0.046/+0.032 SD). The 30–50µm centroid-distance shell is not supported, and the nuclear association is not supported at P14/P21. This is compatible with a very local early growth relationship, not proof of EC-induced CM division. Nuclear/core assignment, lineage uncertainty and late low cycling remain limits.

**The experiment-generating hypothesis is developmental contingency:** vascular apposition may relate more detectably to cell-cycle withdrawal around P7 and to metabolic-regulatory transcription later. This is not a demonstrated timing switch. In particular, low late cycling reduces the ability to detect withdrawal relationships. The selected experiment distinguishes oxygen support, contact and initial CM state rather than assuming that one signal controls a single maturation axis.

## What was reconstructed

| Item | Recovered or implemented |
|---|---|
| Figure 4d population | Exact public `v.CM` labels: 6,298 / 11,086 / 22,532 / 24,588 at P0/P7/P14/P21. |
| Figure 4e selection rule | The released script uses multiple CM labels within four named compact/trabecular and left/right ventricular source categories. This is a different population from Figure 4d. The original anatomical assignments are not exported. |
| Independent ventricular cohorts | Conservative morphology selectors intersected with native tissue territories and the released stage-specific CM label sets: 25,404 / 39,565 / 44,553 / 56,904 CMs. Masks, coordinates and config are reviewable. |
| Maturation calculation | Direct 403-target RNA, pooled per-gene standardized 21-gene mean, CM-reference and within-stage calibration sensitivities, native/core counts, and explicitly feasible 3-bin/100-control sensitivity. None is asserted to be the authors' unreleased `Maturation_Features1` implementation. |
| Neighbor graph | Published-profile symmetric6NN and tissue-screened Delaunay graphs. Rings≥10 are pooled and disconnected profiles flagged. Original graph construction settings/CSV are unavailable. |
| Figure 4e plotting rule | Mean within each distance level followed by the released default R `stats::loess` rule. Both bin counts and individual-cell associations are shown. |
| Physiological outcomes | No perfusion, force, respiration or calcium-function values are imputed. Transcript indices are kept separate from those functions. |

The main text describes averaged standardized expression, while the Methods describe AddModuleScore with 100 matched controls/gene. The control pool/binning/assay are not sufficiently specified. On 403 measured targets, default 24 bins cannot supply 100 distinct controls/bin; a larger assay or different parameters could. Our feasible 3-bin sensitivity is explicitly named and is not proof of an author implementation error or an exact AddModuleScore reproduction.

## Anatomy, sources and interfaces

The independent masks use broad manually specified ventricular morphology selectors, original native cell-territory geometry, local porosity, distance to large voids and outer surface, and 200µm spatial regions. They use no CM maturation values. Dense and cavity-associated classes are **compact-compatible and porous morphology categories**, not recovered LV/RV identities or independently verified trabecular anatomy. Large voids can include section artifacts or vessels; they are not automatically called ventricular lumens. The anatomical masks and limitations are shown directly.

Vascular abundance is measured as nearby capillary count/area within a tissue-screened 50µm centroid neighborhood. Boundary apposition samples 64 equally spaced positions on the original CM polygon and records the fraction within 1/2/5µm of capillary polygons. An outcome-blinded 256-CM-per-age check using 128/256 samples gives high agreement (Spearman generally >0.99 and mean absolute fraction differences <0.005). Refined-boundary and native≤1µm contact counts provide alternative geometries. These operational 2D measurements do not prove membrane protein contact or perfusion.

Source/focal role conflicts are explicit. For example, some independently selected ventricular CMs have public capillary/High.Myo labels. They must not be treated as contacting their own source polygon or as valid zero-distance CM-to-EC examples. Conflicting profiles remain in the atlas and are flagged/excluded from the affected source-specific comparison. This does not establish which identity system is correct. The paper-label CM/source comparison is disjoint by definition.

CM profiles without measured nuclear area remain in the main atlas and eligible main analyses. Nuclear-depth requirements apply only to the separately named nuclear cycling sensitivity.

## Different components, not a single functional maturity scale

The measured 21-gene panel is retained. Additional explicit within-panel groupings distinguish metabolic-regulatory, sarcomere, excitation/contraction and cycling transcripts. A higher “withdrawal” score means lower cycling RNA; it does not establish completed cell-cycle arrest or productive division. The metabolic grouping contains regulatory genes and is not mitochondrial flux.

Broader RNA-only respiratory-chain, sarcomere and calcium-handling panels are also supplied. They do not all increase monotonically. This reinforces the need to distinguish regulation, transcript abundance and actual function. RNA-only genes are not inserted into the Xenium count matrix or presented as measured spatial expression.

Joint models report unadjusted-context, within-region, capture-adjusted and CM-state-adjusted associations separately. Primary associations do not silently condition away receiver expression. Secondary receiver/matrix interaction models are exposed; they do not consistently support a positive matrix gate. Within-section 400µm block intervals, collinearity diagnostics and BH-adjusted approximate block tests accompany coefficients. They do not replace independent hearts. Fixed Ridge comparisons use buffered spatial holdout and show modest component-specific predictive increments.

## RNA ordering and P3/P5/P10

The observed ventricular RNA cohort comprises 43,249 QC nuclei. A 1500-HVG/30-PC neighbor analysis with maturation, validation and cycling genes withheld gives a connected P0-rooted ordering. The pooled age association is approximately 0.90 and disjoint RNA replicate-series associations are approximately 0.89. Each series uses a different RNA library at every age; it is not longitudinal tracking of an animal.

Early gaps remain sparsely bridged: only about 0.6% of P0/P7 neighbor edges cross between those ages, and P7/P14 connections are also sparse. Thus, the ordering is stronger evidence for broad progression than for the path inside a missing interval. Depth remains associated with ordering.

The P3/P5/P10 table contains **conditional transcript-program expectations** from linear and shape-preserving interpolation of observed library endpoints. Intervals combine endpoint-library bootstrap and model choice. They do not measure unknown curvature. The calendar estimates are not a calibrated pseudotime clock, new nuclei, reconstructed spatial sections, missing cell-state distributions or biological replication. Cycling-module means also need not equal the fraction of cycling cells.

## Figures and tables

- [Operational anatomy](01_independent_anatomy_masks.png), [P7 descriptive reconstruction](02_P7_reconstructed_descriptive_relationship.png), [all-stage component/apposition maps](03_all_stage_component_and_apposition_maps.png).
- [Shared anatomy and strict matched comparisons](04_shared_anatomy_and_strict_matched_CM_results.png), [geometry/identity/nuclear sensitivities](05_geometry_identity_and_nuclear_sensitivity.png), [joint-model coefficient gallery](06_joint_context_coefficients.png).
- [Observed RNA ordering](07_observed_RNA_ordering_and_gap_support.png), [conditional missing-day estimates](08_conditional_P3_P5_P10_program_estimates.png).
- [Reconstruction audit](author_reconstruction_audit.tsv), [cohorts](cohort_reconstruction_counts.tsv), [source-role conflicts](focal_source_identity_conflicts.tsv), [outcome definitions](outcome_program_definitions.tsv).
- [Distance-level means](paper_style_neighbor_level_means.tsv), [R LOESS curves](paper_style_R_LOESS_curves.tsv), [descriptive associations](descriptive_proximity_associations.tsv), [covariance partition](shared_anatomy_covariance_partition.tsv).
- [Joint coefficients](anatomy_conditioned_coefficients.tsv.gz), [matched pairs](matched_apposition_CM_pairs.tsv.gz), [all matched outcomes](matched_apposition_component_differences.tsv), [balance checks](matched_interface_covariate_balance.tsv).
- [Contact-definition models](paired_contact_definition_coefficients.tsv.gz), [surface-resolution checks](interface_geometry_sensitivity.tsv), [nuclear cycling coefficients](nuclear_cycling_specificity_coefficients.tsv), [nuclear sensitivity coverage](nuclear_cycling_specificity_coverage.tsv).
- [Spatial holdout](spatially_heldout_mechanism_prediction.tsv), [fold separation checks](spatial_prediction_fold_checks.tsv), [observed RNA library profiles](RNA_observed_library_programs.tsv), [gene support](RNA_program_gene_support.tsv).
- [Observed RNA ordering assessment](RNA_ordering_assessment.tsv), [per-nucleus orderings](RNA_cell_developmental_ordering.tsv.gz), [cross-age graph support](RNA_cross_age_graph_edges.tsv), [P3/P5/P10 expectations](inferred_P3_P5_P10_programs.tsv).
- [Experiment priority and falsifiers](experiment_priority_and_falsifiers.tsv), [structured experiment plan](selected_experiment.json), [numerical validation](numerical_validation.json), [release validation](release_validation.json).

Each stage folder contains original-ID observations, morphology masks, CM component scores, neighbor graphs/levels, interface/matrix context and sensitivity data. Join using `(stage, cell_id)`; cell IDs alone are not globally unique. Original count files and coordinates are unchanged.

## Reproduce

```bash
OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 \
 /workspace/heart-atlas/venv/bin/python analysis/run_cardiac_reconstruction.py \
 --Rscript /workspace/heart-atlas/r-env/bin/Rscript
```

Prerequisites are the existing validated QC/count/segmentation caches and pinned Python dependencies in `analysis/requirements-cardiac-reconstruction.txt`; R 4.5.3 uses base `stats::loess`. Use `--resume` to verify source/code/configuration/artifact hashes and reuse the release. `--steps` permits explicit partial execution with existing prerequisites. No Wormhole retraining is required for this analysis.

Validation independently recomputes measured scores, brute-force source-polygon distances, graph BFS levels, a conditional coefficient, covariance decomposition, strict matching calipers/balance, disjoint RNA series, withheld genes and interpolations. Original count-file hashes are checked against the earlier release. Document checks verify PDF/PNG integrity, links, hashes and file-size limits. Numerical checks do not establish anatomical truth or causality.

Sources: [Wang et al. atlas](https://doi.org/10.1038/s44161-026-00861-z), [released code snapshot](https://github.com/whfbio/PIP_seq/tree/102d80e6ea5de44f594461c727a4b3c6e1bc3508), public WebAtlas labels, GSE290989 and the deposited Xenium datasets; [Lu et al. perinatal vascular study](https://doi.org/10.1038/s41556-021-00747-1) provides the expansion context. Public Xenium derivatives retain attribution to the deposit/article authors, CC BY 4.0. The uploaded journal PDF and original full snRNA expression object are not redistributed here.
