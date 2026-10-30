# Integrated cardiac spatial pipeline

This release establishes a reproducible analysis of **all 408,054 deposited spatial profiles**, their broad/finer identities, segmentation-aware neighborhoods, computational tissue contexts and associations between non-CMs and the measured ventricular CM maturation panel. It identifies candidates for biological follow-up; it cannot establish causal influence from four cross-sectional sections.

[Atlas/analysis PDF](cardiac_spatial_pipeline_report.pdf) · [Configuration](configuration.json) · [Validation](validation.json) · [Model candidates](candidate_relationships.tsv) · [Full joint-model results](joint_association_coefficients.tsv.gz) · [Contact enrichment](contact_enrichment_local_null.tsv) · [Niche stability](niche_model_selection.tsv)

## What is integrated

- The existing reviewed broad atlas provides cardiomyocytes, endothelial cells, fibroblasts, pericytes, smooth-muscle cells, immune cells and epicardial cells. Mixed, ambiguous and QC-flagged profiles remain in the dataset. No source counts, physical coordinates or original broad labels are changed. This is not a complete census of every biological cell type: the panel cannot establish an exhaustive neural/glial population, for example.
- Prioritized CM candidates, reviewed fibroblast families and five endothelial types are retained at their supported resolution. FB/EC primary fine maps require existing strict cell+nuclear support; unresolved fine states remain unresolved. Native RNA-defined substates are not automatically painted onto space.
- Every lineage receives an additional **boundary-core versus nuclear marker audit**. Weak own-lineage support and multiple strong lineage panels are review flags rather than proven doublets. A conservative source/focal-cell mask is tested separately; it never silently replaces the primary labels. Missing nuclear area alone does not remove a CM.
- All primary FB/EC functional exposures now use the same nuclear-anchored count source: FB programs use boundary-core RNA, EC cycling/guidance use nuclear RNA. This replaces the earlier mismatch between native FB model inputs and newer segmentation-supported maps. Whole-cell/native CM scores remain a sensitivity outcome.

## Spatial geometry and neighborhood discovery

**Shared-boundary graph:** profile pairs must share at least 1.275µm of boundary in the existing 0.425µm segmentation grid. Contacts are reconstructed from pixel interfaces, not inferred from expression or centroid proximity. The source segmentation is antibody-guided, nuclear-anchored and retains vendor profile IDs. It does not independently discover missing cells or solve every merge/split error, so a shared boundary is operational adjacency rather than verified biological contact.

**Tissue-screened proximity graphs:** 30/50/100µm centroid pairs are retained only when at least 85% of 21 sampled points along their connecting segment fall within the refined cell-territory union. This reduces straight links across tissue gaps/lumens. It is a sampled barrier screen, not a manually annotated tissue mask or a geodesic-distance solution. At50µm, this rejects approximately2.2%,1.7%,6.2% and9.4% of distance-only links at P0/P7/P14/P21.

Contact enrichment is compared with 199 label permutations **within 200µm spatial blocks**, holding reviewed/unresolved positions fixed. This preserves coarse local composition and goes beyond comparing against each whole section's global cell fractions. The exported p/q values describe this within-section randomization only; they are not animal-level evidence of communication.

Computational contexts use one- and two-step contact-neighbor **resolved lineage proportions**. Focal self-return in two-step walks is removed. Focal RNA, CM maturation, age and global XY are absent from the features. Review burden is held out as a QC diagnostic; training/prediction require at least 3 resolved direct neighbors and≥50% resolved mass in both steps. Training is balanced across stages and excludes reserved 400µm blocks. Gaussian-mixture models test 2–12 components with three seeds/three initializations. Primary selection requires **minimum pairwise seed ARI≥0.8**, then compares BIC among qualifying models. Higher-resolution partitions that fail stability remain evidence of a limitation rather than being promoted to new tissue compartments.

Cluster posterior and seed agreement are computational diagnostics, not independently validated anatomical accuracy. Cell-type distributions and continuous neighborhoods remain usable even where a hard computational context is unresolved.

## Relating non-CMs to ventricular maturation

The inherited strict ventricular cohort contains 149,623 profiles, selected partly by measured Myl2. This selection is not independent of the maturation panel; omitting Myl2 from an outcome does not recover excluded low-Myl2 CMs. This limitation remains explicit.

The primary outcome is the measured **21-gene mean log1p(panel CP10k) from boundary-core RNA**. The denominator is the 403-gene assay, not the whole transcriptome. Sensitivities include native-profile 21-gene scores, a 20-gene score omitting Myl2, the four supplied maturation categories, and an equal-category standardized score. The latter is relative within age and should not be interpreted as an absolute cross-age scale. All outcome definitions are frozen before ranking spatial effects.

Primary joint models use 50µm tissue-screened neighborhoods with≥5 resolved neighbors,≥1 eligible EC and≥1 eligible FB, adequate core CM RNA excluding panel outcomes, and available anatomical proxies. This yields 13,588/31,830/30,300/29,031 CMs at P0/P7/P14/P21. **The same CMs and joint exposures are used across adjustment levels within each geometry**. Contact analyses use fewer CMs; matched contact/50µm comparisons explicitly restrict both to the identical intersection, rather than attributing cohort differences to geometry. Coverage and missing-nucleus retention are exported.

Six non-CM composition fractions and ten prespecified measured programs are tested jointly. These include FB fibrillar/basement/elastic/remodeling RNA, EC cycling and guidance, and limited immune/mural/pericyte/epicardial proxies. [Actual genes and interpretation limits](program_definitions.tsv) are exported. No maturation-panel gene is included in these exposures. A proxy list is not automatically a validated coexpression program; contractility/vascular-support proxies overlap lineage markers. Strong collinearity and insufficient exposure variation are flagged.

The three adjustment levels answer progressively narrower questions:

1. **Technical:** CM capture excluding score genes, profile area, nuclear-area missingness, neighborhood density/review fraction, source coverage/capture and source-presence indicators.
2. **Anatomical:** adds spatial splines, a spatial interaction, and existing arterial-envelope/endocardial-proxy distances.
3. **CM state:** also conditions on CM cycling and Slit2 expression. This can remove a shared developmental pathway and should not be read as the total biological effect of non-CMs.

Composition coefficients represent relative changes in the neighborhood with CMs as the omitted composition component; they are not independent changes in absolute cell numbers. Measured RNA programs do not measure deposited matrix, mechanics, secreted protein, receptor activation or earlier exposure.

Primary intervals use 512 resamples of 400µm spatial blocks; sensitivity models use 128. Buffered five-fold prediction uses the same 400µm blocks and removes training CMs within 150µm of test CMs—more than twice the primary exposure radius. Baseline/composition/program/joint Ridge models and a fixed-hyperparameter nonlinear boosted-tree sensitivity are compared on held-out core21-gene scores. A small added R² means little extra prediction beyond that baseline, not biological irrelevance. Bootstrap intervals condition on the current identities/count views and do not propagate every reference or segmentation uncertainty. One section per age cannot provide biological-replicate inference.

## Findings worth pursuing

- Smooth-muscle, epicardial and immune profiles show spatial concentration, including homotypic contact enrichment after preserving coarse local composition. Some heterotypic relationships, such as early fibroblast/immune contacts, are section/stage-specific hypotheses. [Complete local-null table](contact_enrichment_local_null.tsv) prevents selection of only favorable patterns.
- EC cycling and the measured guidance proxy do not show a consistent negative relationship with the overall CM maturation score at P0/P7 in the harmonized joint models. This does not establish absence of an angiogenic influence, and a sprouting-specific phenotype remains incompletely measured.
- Higher local fibroblast proportion has a **small negative association** with overall CM maturation at later stages. P21 meets the operational robustness screen, including same-direction conservative-QC and paired-geometry checks. This is a relative-composition association, not proof that fibroblasts suppress maturation; individual matrix programs provide weaker conditional evidence. Full component and geometry results remain available.
- Candidate signaling tables provide a **curated RNA-expression hypothesis screen**, not CellChat output or inferred interaction probabilities. They use 13 ligand/receptor-component pairs, per-library corrected shared-panel pseudobulks,≥10 source/receiver nuclei and expression support in both libraries. Receptor complexes/co-receptors, protein signaling and functional effects are not established. The threshold uses shared 381-target CPM, not whole-transcriptome CPM. No signaling mechanism is validated by expression alone.

The candidate screen requires a primary effect≥0.05SD, interval direction, native/category-balanced/Myl2-withheld score sign agreement, paired contact/proximity sign agreement, technical-model sign agreement, adequate partial exposure variance and conservative-QC support. The0.05 threshold is an operational ranking choice, not a validated biological effect threshold. Multiple models are explored; these gates do not constitute multiplicity-corrected discoveries or independent replication.

## Do we need CellCharter or Wasserstein Wormhole?

**Neither is required for this primary pipeline.** CellCharter is a reasonable complementary niche-discovery benchmark. The actual default here is an interpretable scikit-learn GMM on graph-derived compositions; it is not a CellCharter analysis. An explicit optional adapter runs CellCharter's clustering backend on the same exported features. A full transcriptomic CellCharter comparison would require a separately assessed expression representation and would answer a broader state-discovery question.

The uploaded **Haviv et al., Wasserstein Wormhole: Scalable Optimal Transport Distance with Transformers, ICML 2024/arXiv 2404.09411** embeds entire neighborhood distributions so Euclidean distances approximate OT distances. Section 6.5 is directly relevant to cellular niches. It is not cell segmentation, snRNA–Xenium allocation, or causal inference. We export observed non-CM type/program point clouds, benchmark 6,000 neighborhood pairs with sliced OT and 256 small exact OT solves, and provide an optional verified 0.3.7 training adapter. The pilot tests distance agreement only and uses restricted module/type features; it does not test the full 403-gene potential of Wormhole or show improved CM prediction. The transformer was **not trained** in this release.

The current environment does not install the optional PyTorch/CellCharter or JAX/Wormhole stacks. Separate environments preserve the primary analysis dependencies. Official sources: [CellCharter](https://github.com/CSOgroup/cellcharter), [Wormhole paper](https://arxiv.org/abs/2404.09411), [official package/API](https://pypi.org/project/wassersteinwormhole/0.3.7/). Method versions/API checks and the uploaded paper hash are in [method assessment](method_assessment.json).

## Outputs and reproduction

| Stage | Integrated atlas | Supported subtypes | All-profile annotations | Analysis-ready niche object |
|---|---|---|---|---|
| P0 | [PNG](P0_integrated_atlas.png) | [PNG](P0_supported_subtypes.png) | [TSV.gz](P0/atlas_annotations.tsv.gz) | [H5AD](P0/niche_analysis_ready.h5ad) |
| P7 | [PNG](P7_integrated_atlas.png) | [PNG](P7_supported_subtypes.png) | [TSV.gz](P7/atlas_annotations.tsv.gz) | [H5AD](P7/niche_analysis_ready.h5ad) |
| P14 | [PNG](P14_integrated_atlas.png) | [PNG](P14_supported_subtypes.png) | [TSV.gz](P14/atlas_annotations.tsv.gz) | [H5AD](P14/niche_analysis_ready.h5ad) |
| P21 | [PNG](P21_integrated_atlas.png) | [PNG](P21_supported_subtypes.png) | [TSV.gz](P21/atlas_annotations.tsv.gz) | [H5AD](P21/niche_analysis_ready.h5ad) |

Use `stage_cell_id` across stages. The niche H5AD has annotation/coordinates/features/contact graph and deliberately **no X RNA matrix**; original observed count files are referenced through source paths/hashes. Per-stage sparse NPZ graphs, contact-edge lengths, barrier-audit samples, neighborhood exposures, QC flags, OT point clouds/subsamples and ROI provenance are included. [Atlas counts](atlas_counts.tsv), [geometry diagnostics](graph_summary.tsv), [cohort coverage](model_cohort_coverage.tsv), [QC sensitivity](conservative_QC_association_sensitivity.tsv.gz), [prediction](spatial_prediction.tsv), [nonlinear prediction](nonlinear_prediction.tsv), [signaling hypotheses](candidate_signaling_expression.tsv), [transport pilot](transport_pilot_summary.tsv).

```bash
OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 \
  /workspace/heart-atlas/venv/bin/python analysis/run_cardiac_spatial_pipeline.py
# Reuse only if configuration, inputs, code and artifacts still match:
/workspace/heart-atlas/venv/bin/python analysis/run_cardiac_spatial_pipeline.py --resume
# Optional backends, in separate environments; not run by the default command:
python analysis/cardiac_spatial_niches.py --step cellcharter
python analysis/cardiac_spatial_niches.py --step wormhole
```

Copy the JSON configuration with a different `output` directory for an isolated alternative run; pass it with `--config`. Explicit `--steps` supports partial execution after prerequisites exist. Runtime versions, source hashes and artifact/script hashes support review. System dependency: Poppler `pdfinfo` for PDF validation. Validation independently reconstructs sampled complete contact boundaries, tissue screens, observed scores, neighborhood means, model coefficients, graph-only niche features/posteriors and small exact OT distances. These implementation checks do not validate anatomy or causality.

Source counts, segmentation and deposited Xenium data retain attribution to the paper's authors as documented in the preceding atlas/segmentation analyses. The uploaded methodological paper is cited and hashed, not redistributed.
