# Persistent neighborhoods and changing stromal–vascular states in postnatal heart development

[Editable Word manuscript](Postnatal_cardiac_development_manuscript.docx) · [Rendered manuscript PDF](Postnatal_cardiac_development_manuscript.pdf) · [Six integrated figure plates](Cardiac_development_integrated_figures.pdf) · [Publication figure bundle](Cardiac_development_publication_figures.zip) · [Source-data bundle](Cardiac_development_figure_source_data.zip)

**Central message:** postnatal cardiac development is consistent with remodeling of recurring regional tissue neighborhoods. Fibroblast matrix specialization, vascular growth and CM developmental programs change together, but need not follow a single synchronized maturation axis.

This is a manuscript-style **secondary-analysis draft**, with abstract, introduction, evidence-linked results, discussion, methods, references, six multi-panel figures, full legends and an experiment roadmap. It distinguishes observed associations, candidate subtype continuity and proposed mechanisms. It does not present the recovered public labels as newly discovered types or claim independent biological replication.

## The figure sequence

| Figure | Biological message | Files |
|---|---|---|
| 1 | Recurring myocardial, vascular and stromal/surface organization anchors comparison | [PNG](figures/Figure_1_persistent_spatial_architecture.png) · [PDF](figures/Figure_1_persistent_spatial_architecture.pdf) · [SVG](figures/Figure_1_persistent_spatial_architecture.svg) |
| 2 | Persistent regional CM identity overlaps with partial peri-arterial response programs | [PNG](figures/Figure_2_regional_cardiomyocyte_programs.png) · [PDF](figures/Figure_2_regional_cardiomyocyte_programs.pdf) · [SVG](figures/Figure_2_regional_cardiomyocyte_programs.svg) |
| 3 | Fibroblast programs specialize within overlapping interstitial environments; RNA resemblance nominates continuity | [PNG](figures/Figure_3_fibroblast_state_remodeling.png) · [PDF](figures/Figure_3_fibroblast_state_remodeling.pdf) · [SVG](figures/Figure_3_fibroblast_state_remodeling.svg) |
| 4 | EC growth and capillary apposition relate to different CM developmental components | [PNG](figures/Figure_4_vascular_growth_and_CM_components.png) · [PDF](figures/Figure_4_vascular_growth_and_CM_components.pdf) · [SVG](figures/Figure_4_vascular_growth_and_CM_components.svg) |
| 5 | Transport compares remodeling neighborhoods but does not establish a unique fate map | [PNG](figures/Figure_5_transport_constrained_correspondence.png) · [PDF](figures/Figure_5_transport_constrained_correspondence.pdf) · [SVG](figures/Figure_5_transport_constrained_correspondence.svg) |
| 6 | A regional tissue-unit model yields oxygen/contact/CM-state and lineage-tracing experiments | [PNG](figures/Figure_6_developmental_model_and_tests.png) · [PDF](figures/Figure_6_developmental_model_and_tests.pdf) · [SVG](figures/Figure_6_developmental_model_and_tests.svg) |

All plots were assembled anew from deposited data and frozen reanalysis tables. PNGs are 300 dpi; PDFs preserve vector axes/text, with dense spatial point layers rasterized, and SVGs retain editable text/graphics. Maps preserve original spatial coordinates and use physical scale bars. The Word file embeds the full current figure plates; detailed legends follow the plates. A layout-rendered PDF is provided for review.

## Additional analyses performed for this story

- **Reviewed fibroblast matrix and marker profiles:** library×family/substate pseudobulks from the current clean RNA reference, using full-RNA denominators and both original and ambient-corrected counts. The late Col8a1-versus-Gfpt2 basement-matrix difference is supported in both libraries and count views. P21_F0 retains its reviewed Col8a1/Svep1 identity.
- **Candidate state continuity:** frozen gene-withheld RNA signature matches are relabeled using the current reviewed families. Selected subgroup matches nominate tracing experiments; dashed arrows are not transport probabilities or lineage ancestry.
- **Restricted exact optimal transport:** 96 reserved observed CM-centered neighborhoods, 24 per age and 12 per endocardial-distance category, within the existing CM/EC-enriched context. Exact non-CM cloud costs precede an anatomy-proxy-restricted uniform-mass neighborhood correspondence. All 864 adjacent-stage pair costs and 528 within-section pair costs are retained.
- **Wormhole and correspondence sensitivity:** the saved encoder is reused, not retrained. Broad distance ranking can be reasonable while fine nearest-neighborhood and transport-plan agreement are limited. State/lineage weighting, source count, mean features and entropic versus unregularized correspondence are exposed. Uniform mass is not observed expansion or fate probability.

The transport cost decomposition for composition-similar pairs is conditional on the chosen feature scaling and selection. It is not variance explained, causal importance or proof of a biological transition. Within-section variation overlaps adjacent-stage differences. Residual section/capture/anatomical effects remain possible.

## Source data and validation

- [Fibroblast library matrix programs](fibroblast_library_matrix_programs.tsv), [library marker profiles](fibroblast_library_marker_profiles.tsv), [program genes](fibroblast_RNA_program_definitions.tsv), [candidate continuity table with subgroup IDs](reviewed_fibroblast_candidate_state_continuities.tsv).
- [Selected observed neighborhoods](OT_selected_observed_neighborhoods.tsv), [all adjacent-stage pair costs](anatomy_stratified_neighborhood_pair_costs.tsv), [within-section baseline](within_section_neighborhood_distance_baseline.tsv), [correspondence plans](adjacent_stage_neighborhood_correspondence_plans.tsv), [metric sensitivities](neighborhood_transport_metric_sensitivity.tsv), [conditional cost decomposition](composition_similar_neighborhood_cost_decomposition.tsv), [design](OT_design.json).
- [Panel/message key](figure_message_and_panel_key.tsv), [editable manuscript text](Manuscript_text.md), [numerical checks](numerical_validation.json), [document checks](document_validation.json), [release checks](release_validation.json).

Numerical checks independently reconstruct six cloud LP costs with dense constraints and six neighborhood LP objectives with Hungarian assignment; verify all marginals, decompositions, source membership, selected-cohort spacing/eligibility, withheld outcomes and metric summaries; and directly recompute two RNA pseudobulks. Document validation checks that the Word file embeds the exact current figure PNGs, renders six full-page figure plates within page bounds and includes valid 300-dpi PNG, vector PDF and editable SVG exports. Source/code/artifact hashes support verified reuse. These checks establish computational consistency, not anatomical truth, ancestry or causality.

## Reproduce

Prerequisites are the existing QC, segmentation and count caches referenced by upstream manifests. Use the scientific environment specified by `analysis/requirements-cardiac-reconstruction.txt`; the optional document environment is specified by `analysis/requirements-cardiac-story-document.txt`. LibreOffice/soffice supplies Word-to-PDF rendering.

```bash
OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 \
 /workspace/heart-atlas/venv/bin/python analysis/run_cardiac_story.py \
 --document-python python

# Verify inputs, scripts and artifacts before reusing the completed release:
/workspace/heart-atlas/venv/bin/python analysis/run_cardiac_story.py --resume
```

`--steps` permits explicit partial execution; dependencies must already exist. Narrative and interpretation are deliberately authored in `analysis/write_cardiac_story.py` and must be reviewed if analysis definitions change. Intermediate P3/P5/P10 predictions remain in the upstream reconstruction and are not shown as observed data here.

## Scope and attribution

There is one available spatial section per age. RNA libraries are sensitivity units and are not assumed independent hearts. Some fine cell identities remain unresolved, particularly sparse fibroblast/EC subtypes; immune, mural and epicardial fine spatial taxonomy is incomplete. No true temporal lineage, absolute expansion, oxygen/perfusion, matrix mechanics or CM functional outcome is inferred as directly measured. Original author masks, score-bearing objects and neighbor inputs remain unavailable for exact Figure 4e replication.

Data: [GSE290989](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE290989), Xenium [P0/P7](https://zenodo.org/records/14961929) and [P14/P21](https://zenodo.org/records/14961942), and [Wang et al.](https://doi.org/10.1038/s44161-026-00861-z). Biological context: [Lu et al.](https://doi.org/10.1038/s41556-021-00747-1). Transport method: [Haviv et al.](https://arxiv.org/abs/2404.09411). Public data derivatives retain original article/deposit attribution (CC BY 4.0). No uploaded journal PDF is redistributed.
