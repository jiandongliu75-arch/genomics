# Subtype-first cardiac neighborhoods and within-stage CM maturation

[Revised Word manuscript](Subtype_resolved_cardiac_development_manuscript.docx) · [Rendered manuscript](Subtype_resolved_cardiac_development_manuscript.pdf) · [Six integrated figure plates](Subtype_resolved_integrated_figures.pdf) · [Publication figures](Subtype_resolved_publication_figures.zip) · [Source tables](Subtype_resolved_figure_source_tables.zip)

**The primary unit is stage × CM identity × within-CM-state maturation × specific fibroblast/endothelial source subtypes.** The earlier two broad lineage-composition contexts are background descriptions, not a limit on known subtype analyses. No unsupervised niche-cluster stability filter is imposed here.

## What changed

- **CM identity and maturation are kept separate.** Ankrd1, working ventricular, Slit2, Ddc and cycling CM states are retained where supported by the ventricular gate. Continuous measured scores and within-stage×CM-state quartiles represent different maturation levels inside each identity, including P0/P7.
- **Neighbors are specific FB/EC types.** Postn, Gfpt2, Col8a1 and other supported fibroblast families are modeled alongside capillary, endocardial, arterial and other supported EC types. Unknown subtype sources remain in lineage denominators. Sparse source calls are not biological absence.
- **Regional and residual local relationships are distinguished.** Technical/composition, anatomical-proxy and within-200µm-region models answer different questions. Composite and metabolic/sarcomere/excitation-contraction/cycling-withdrawal RNA components are retained.
- **CM-state specificity is tested.** Ankrd1-versus-working slope differences use shared spatial-block bootstrap weights and common units rather than comparing separate p-values.
- **Geometry is checked on the same CMs.** New 50/30µm and 50µm/contact intersections separate geometry sensitivity from focal-cohort selection.
- **Transport respects CM identity.** Ninety-six reserved neighborhoods are compared within the same CM state and endocardial-distance category. Known source subtypes and measured CM maturity remain explicit descriptors. There is no two-context cluster filter and no claim that the old encoder was retrained on new subtype labels.

## Current results and limits

Ankrd1 and working CMs each span substantial maturation-score distributions at P0/P7. Raw source-subtype composition differs by CM identity and maturation group, but many overall-score associations in these two CM states become weak after local controls. This distinguishes regional organization from a residual local mechanism.

Two more specific candidates merit follow-up: Postn-supported FB context accompanies lower overall/metabolic-regulatory RNA within P7 cycling CMs; endocardial-like EC context accompanies higher sarcomere RNA within P14 Ankrd1-CMs. These are exploratory within-section findings. Stricter source definitions preserve some directions, broader public annotations weaken certainty, and the P7 result is not supported at30µm. Contact tests are smaller and do not establish direct contact signaling. The full result grid is retained.

Source coverage is strongly family/stage dependent. P0/P21 have only one well-supported primary FB family under current thresholds; comparisons there partly concern supported-versus-untyped FB context, not clean competition among all RNA families. Unknown labels are not filled with public or soft proposals. Public source labels are a separately named comparator.

## Integrated figures

| Figure | Message | Exports |
|---|---|---|
| 1 | CM identity and maturation are distinct within P0/P7 | [PNG](figures/Figure_1_CM_identity_and_maturation.png) · [PDF](figures/Figure_1_CM_identity_and_maturation.pdf) · [SVG](figures/Figure_1_CM_identity_and_maturation.svg) |
| 2 | Specific FB/EC subtype profiles around each CM identity and maturity group | [PNG](figures/Figure_2_known_subtype_neighborhoods.png) · [PDF](figures/Figure_2_known_subtype_neighborhoods.pdf) · [SVG](figures/Figure_2_known_subtype_neighborhoods.svg) |
| 3 | Component- and CM-state-specific conditional associations and sensitivities | [PNG](figures/Figure_3_CM_state_specific_associations.png) · [PDF](figures/Figure_3_CM_state_specific_associations.pdf) · [SVG](figures/Figure_3_CM_state_specific_associations.svg) |
| 4 | Reviewed fibroblast specialization and possible RNA-state continuity | [PNG](figures/Figure_4_fibroblast_specialization.png) · [PDF](figures/Figure_4_fibroblast_specialization.pdf) · [SVG](figures/Figure_4_fibroblast_specialization.svg) |
| 5 | CM-identity-conditioned transport and held-out maturity audit | [PNG](figures/Figure_5_CM_identity_conditioned_transport.png) · [PDF](figures/Figure_5_CM_identity_conditioned_transport.pdf) · [SVG](figures/Figure_5_CM_identity_conditioned_transport.svg) |
| 6 | Subtype-defined developmental hypotheses and matched experiments | [PNG](figures/Figure_6_subtype_developmental_framework.png) · [PDF](figures/Figure_6_subtype_developmental_framework.pdf) · [SVG](figures/Figure_6_subtype_developmental_framework.svg) |

All new plots use original physical coordinates. Figures4 reuses the explicitly versioned earlier fibroblast plate; Figures1–3/5–6 are new. PNGs and rasterized PDF layers are300dpi; text/axes remain vector in PDF and editable in SVG. Full legends are embedded in the Word manuscript. Source masks and count views are stated rather than silently combined.

## Main tables and checks

- [CM cohorts](CM_state_cohort_coverage.tsv), [source subtype coverage](source_subtype_coverage.tsv), [source fractions by CM identity×maturation](CM_state_maturation_subtype_neighborhood_means.tsv), [P0/P7 score distributions](P0_P7_CM_identity_maturation_distributions.tsv).
- [All joint coefficients](subtype_conditioned_maturation_coefficients.tsv), [CM-state slope differences](Ankrd1_versus_working_CM_slope_differences.tsv), [same-CM geometry](same_CM_geometry_sensitivity.tsv), [high/low matched pairs](within_CM_state_maturation_matched_pairs.tsv.gz), [matched subtype differences](matched_maturation_subtype_neighborhood_differences.tsv).
- [CM-identity-restricted OT neighborhoods](CM_identity_restricted_OT_neighborhoods.tsv), [all pair costs](CM_identity_restricted_OT_pair_costs.tsv), [correspondence plans](CM_identity_restricted_OT_correspondence.tsv), [Wormhole assessment](CM_identity_restricted_Wormhole_assessment.tsv), [held-out maturity audit](CM_maturation_rank_correspondence_audit.tsv), [within-section distances](CM_state_within_section_neighborhood_distances.tsv).
- [Numerical validation](numerical_validation.json), [document validation](document_validation.json), [release checks](release_validation.json), [readable manuscript text](Manuscript_text.md).

Each stage has all eligible CM neighborhood records for proximity50, proximity30 and contact graphs, preserving native IDs/coordinates, source-count denominators, all subtype views, coverage, score ranks and cohort flags. Join by stage+cell ID. These are source-derived annotations; original403-target matrices are unchanged.

## Reproduce

Use the existing scientific environment and frozen upstream caches. Document dependencies are listed in `analysis/requirements-cardiac-story-document.txt`; LibreOffice/soffice renders the Word file.

```bash
OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 \
 /workspace/heart-atlas/venv/bin/python analysis/run_cardiac_subtype_story.py \
 --document-python python

/workspace/heart-atlas/venv/bin/python analysis/run_cardiac_subtype_story.py --resume
```

Validation independently reconstructs neighbor subtype fractions and score ranks, checks disjoint source/focal roles, refits target coefficients by a separate centered least-squares path, verifies matched calipers/cohort intersections and CM-state slope arithmetic, reconstructs exact OT costs with dense constraints and outer objectives by Hungarian assignment, checks marginals and confirms the Word file embeds the exact current figure plates. These are implementation checks rather than anatomical, lineage or causal validation.

One spatial section per age remains the main replication limit. RNA libraries are sensitivity units, not assumed independent hearts. CM state labels here explicitly use public final spatial identities intersected with our ventricular/QC gates. Known RNA subtypes can exceed confident spatial resolution. The saved Wormhole model is not a new fine-label-trained encoder, and transport mass is not a fate probability, expansion estimate or temporal rate. Candidate subtype relationships and transitions require biological replication and functional/lineage experiments.

Data attribution: [Wang et al.](https://doi.org/10.1038/s44161-026-00861-z), [GSE290989](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE290989), Xenium Zenodo [P0/P7](https://zenodo.org/records/14961929) and [P14/P21](https://zenodo.org/records/14961942), public WebAtlas labels; public derivatives retain article/deposit attribution(CC BY4.0). Context and methods: [Lu et al.](https://doi.org/10.1038/s41556-021-00747-1), [Wasserstein Wormhole](https://arxiv.org/abs/2404.09411). Uploaded journal PDFs are not redistributed.
