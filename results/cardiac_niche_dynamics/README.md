# Cardiac niches, vascular growth and CM maturation

This extension reviews the two papers against the existing QC atlas, tests measured vascular states and segmentation sensitivity, and actually trains and evaluates Wasserstein Wormhole. The earlier `cardiac_spatial_pipeline` release remains intact as the comparison.

**Report:** [13-page PDF](cardiac_niche_coordination_report.pdf). Figures are also individual PNG files for browser viewing.

## Findings

- The two seed-stable composition contexts describe **different neighbor mixtures**: fibroblast/mural enriched versus CM/EC enriched. They are not two new cell types or a claim that the heart contains only two niches. The authors used expression-derived scVI features and CellCharter with 12 clusters; our previous composition-only analysis used different inputs.
- Smooth-muscle, epicardial and immune self-concentration persists with native-polygon and refined-boundary graphs. P0 fibroblast–immune contacts are approximately **1.63-fold native / 1.72-fold refined** above a local abundance-matched shuffle expectation. This supports local organization, not signaling direction or a causal maturation effect.
- In the corrected whole-transcriptome RNA data, EC cycling is higher at P0/P7 than P14/P21 in the two libraries per age. A six-gene sprouting panel is much less convincing: most genes are sparsely detected, and it does not form a coherent positive module after depth/family adjustment. Its mean cannot be interpreted as an angiogenesis rate.
- Eligible spatial ECs with both **Flt4 and Dll4 core detection**, not assigned endocardial or lymphatic nuclear identity, number **383/3,732 at P0, 1,664/8,944 at P7, 972/8,086 at P14, 1,079/8,601 at P21**. These are measured growth-compatible profiles, not morphologically verified sprouts. Many eligible ECs have unresolved nuclear subtype; nonassignment does not establish coronary origin.
- There is **no consistent negative relationship at both P0 and P7** between EC growth/cycling proxies and the overall CM 21-gene maturation index. Relationships vary with EC definition, RNA view and geometry.
- The initial P7 50µm model associates neighboring EC nuclear **Ptprm** with slightly higher CM maturation (about **+0.015 SD per exposure SD**, receiver-conditioned). On **11,619 identical CMs**, with source capture/coverage and receiver adjustment, the coefficients are **+0.0189 [0.0020, 0.0346]** at 50µm, **+0.0122 [−0.0008, 0.0232]** for refined boundary contact, and **+0.0096 [−0.0042, 0.0207]** for native polygons within 1µm. These are small conditional associations; the boundary comparisons do not establish a contact mechanism.
- In a separate identical-cohort P7 source-label comparison, the EC sender estimate is **+0.0147 [0.0066, 0.0234]** for primary ECs, **+0.0153 [0.0066, 0.0242]** for all broad ECs, and **+0.0028 [−0.0054, 0.0117]** when the sender population is the authors’ capillary labels including `cap.EC.High.Myo.`. This highlights population-definition sensitivity, not evidence that either annotation is ground truth.
- **Wormhole was trained**, but it does not pass the all-age metric gate and adds little spatially held-out maturation prediction. It is not adopted as the primary representation. Exact empirical OT between descriptive mean-state archetype medoids is retained as an exploratory cross-sectional comparison, never a reconstruction of cell lineage or temporal rates.

## Why the two papers change the interpretation

[Lu et al., Nature Cell Biology (2021)](https://doi.org/10.1038/s41556-021-00747-1) traced a pre-existing inner myocardial coronary plexus into neonatal coronary vessels, while late fetal/neonatal endocardium did not measurably contribute in the tested tracing windows. Earlier embryonic endocardial contributions are a separate question. Flt4/VEGFR3 and Dll4 identify the growth-associated plexus, but **Dll4 deletion enhances expansion and NICD activation suppresses it**: DLL4/Notch is feedback restraint, not a simple increasing sprouting score. Regional hypoxia, proliferation and vessel morphology were measured with experiments, not inferred from a few transcripts. This paper did not directly test the CM maturation hypothesis.

[Wang et al., Nature Cardiovascular Research (2026)](https://doi.org/10.1038/s44161-026-00861-z) reports a **positive P7** relationship of capillary/proliferating EC and pericyte proximity to CM maturation. It proposes PTPRM–PTPRM contact signaling, with boundary immunostaining and CM-directed perturbation. Those experiments support a PTPRM role but do not independently establish that the necessary sender is an EC rather than another CM. Vessel growth, later vascular integration, and CM maturation can be coordinated without sprouting always marking locally immature CMs.

The full uploaded 13-page Lu PDF and the full cached PMC13561863 text/XML were reviewed. Figure- and method-specific implications are in [the paper review table](paper_review_and_pipeline_changes.tsv) and [METHODS_REVIEW.md](METHODS_REVIEW.md).

## View the figures

| Figure | Contents |
|---|---|
| [Contexts and concentration](01_contexts_and_concentration.png) | Neighbor mixtures and native/refined local contact enrichment |
| [RNA programs](02_RNA_vascular_programs.png) | Two libraries per age; cycling versus sparse sprouting support |
| [Measured EC states](03_measured_vascular_state_maps.png) | All four ages on original section coordinates |
| [Early CM/EC maps](04_early_CM_vascular_maps.png) | Measured CM maturation, growth-compatible ECs and EC nuclear PTPRM |
| [Graph-ring profiles](05_graph_ring_maturation.png) | Native/core scores and alternative EC definitions |
| [Initial adjusted models](06_conditional_growth_and_contact_associations.png) | Small growth/cycling/PTPRM associations |
| [Same-CM geometry sensitivity](07_paired_geometry_depth_sensitivity.png) | Source capture, coverage, receiver and geometry controls |
| [Source-definition sensitivity](08_source_selection_sensitivity.png) | High myocardial EC audit and alternate sender populations |
| [Wormhole distance validation](09_Wormhole_distance_benchmark.png) | Held-out transport fidelity versus neighborhood averages |
| [Wormhole phenotype prediction](10_Wormhole_spatial_prediction.png) | Buffered spatial holdout on the same CM cohort |
| [Exploratory correspondence](11_cross_sectional_state_correspondence.png) | Exact medoid OT; separate sections, not cell trajectories |

## Tables and reusable data

- [Measured-panel coverage](vascular_program_panel_coverage.tsv), [RNA library programs](RNA_program_library_profiles.tsv), [gene detection](RNA_measured_and_unmeasured_gene_support.tsv), [module coherence](RNA_vascular_module_coherence.tsv).
- [EC definition and distance summaries](vascular_stage_summary.tsv), [graph coverage/censoring](ring_graph_coverage.tsv), [ring correlations](graph_ring_correlations.tsv).
- [Initial coefficients](vascular_maturation_coefficients.tsv.gz), [same-CM geometry coefficients](paired_geometry_coefficients.tsv.gz), [P7 sender-source sensitivity](P7_Ptprm_source_sensitivity.tsv.gz), [paired cohort sizes](paired_geometry_cohort.tsv).
- [Public EC label/source audit](author_EC_state_source_selection_audit.tsv), [nuclear subtype coverage](EC_source_nuclear_subtype_coverage.tsv), [native/refined contact enrichment](native_refined_architecture_enrichment.tsv), [coarse-context compositions](coarse_context_composition_explanation.tsv).
- [Feasible maturation controls](feasible_allcell_matched_controls.tsv), [reported-control feasibility audit](reported_maturation_control_feasibility.tsv).
- [Wormhole execution](Wormhole_execution.json), [distance benchmark](Wormhole_reserved_distance_benchmark.tsv), [metric assessment](Wormhole_metric_assessment.tsv), [spatial prediction](Wormhole_maturation_prediction.tsv), [training loss](Wormhole_training_loss.tsv), saved parameters and held-out/train cloud IDs.
- [Exploratory archetypes](cross_sectional_state_archetypes.tsv), [cell assignments](cross_sectional_state_assignments.tsv.gz), [cross-age medoid costs](cross_age_state_correspondence.tsv), [interpretation gate](dynamic_interpretation_gate.json).
- Each `P0/`, `P7/`, `P14/`, `P21/` directory contains per-profile measured states, CM outcomes/exposures, ring distances and censor flags, native boundary graphs at 0.5/1/2µm, pruned Delaunay graphs, non-CM state points/cloud memberships and learned embeddings. Indices preserve original IDs. Join stages using `(stage, cell_id)`, not the local cell ID alone.

## Wormhole benchmark

| Age | Wormhole / Sinkhorn Spearman | Mean features / Sinkhorn Spearman | Metric gate |
|---|---:|---:|---|
| P0 | 0.819 | 0.902 | Fail |
| P7 | 0.924 | 0.911 | Pass |
| P14 | 0.758 | 0.929 | Fail |
| P21 | 0.801 | 0.921 | Fail |

Training uses 3,200 neighborhood distributions (800/age), 800 reserved assessment distributions, up to 32 non-CMs/cloud, 370 observed non-CM genes reduced to 16 PCs plus six-lineage composition, a 32-dimensional encoder, two transformer blocks/two heads/hidden size64, 1,200 steps, batch8, and Sinkhorn epsilon0.1/100 iterations. This is a CPU benchmark, not a reproduction of the uploaded paper’s default architecture. The package trains **mean squared embedding coordinate distance**, which is used for the distance benchmark. Exact small-cloud W² LP costs provide an independent check.

All CM maturation genes are excluded from the source state features **and their normalization denominator**. PCA excludes every source appearing in reserved test clouds. Training and test CM regions have a ≥150µm buffer, so their 50µm source clouds do not overlap. The maturation outcome never trains the encoder. Fixed Ridge prediction adds at most ~0.0011 R² using Wormhole relative to the baseline and performs worse than simple mean states at P0/P7/P14. This benchmark does not show useful added biological inference from Wormhole.

The six archetypes per age are an explicitly fixed, exploratory summary of mean state features following the failed all-age metric gate. They are **not** seed-stable cell niches, the authors’ 12 niches, tracked cells, or probability-conserving developmental flows.

## Reproduce and validate

Requires the existing validated v1 spatial pipeline and its local count/geometry caches; original data were not duplicated. See upstream source manifests for Zenodo/GEO provenance. Use Python3.12 and pinned primary dependencies in `analysis/requirements-cardiac-niche-dynamics.txt`; use a separate environment for `analysis/requirements-cardiac-wormhole.txt`.

```bash
OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 \
  /workspace/heart-atlas/venv/bin/python analysis/run_cardiac_niche_dynamics.py \
  --wormhole-python /workspace/heart-atlas/wormhole-venv/bin/python
```

To verify a completed release without retraining, use `--resume`. Partial execution is available with `--steps`; prerequisites must already exist. `fit_wormhole_cardiac_niches.py --step benchmark` can assess saved embeddings without repeating training. The runner verifies source/code/artifact hashes before reuse. The default bootstrap helper is inherited from v1 with seed731, 400µm blocks and 256 draws; the new sampling/PCA/Wormhole seed is947.

[Validation](validation.json) checks independent native-polygon distances, BFS graph rings, measured count/score arithmetic, source-capture arithmetic, same-CM cohort sizes, a coefficient by independent residualization, held-out source exclusion, actual training-loss reduction, metric gates, exact OT with independent dense LP constraints, PDF/PNG integrity and provenance hashes. Numerical checks do not establish biological segmentation accuracy.

**Limits:** one spatial section per age; sparse and selected panel; strict source selection and unresolved EC families; residual segmentation/capture effects; approximate anatomy; ventricular selection partly uses Myl2, one of the maturation-panel genes. The earlier release’s Myl2-withheld/category-balanced scores remain relevant sensitivity analyses. No local oxygen/perfusion measurement, lineage tracing, EC-specific perturbation or independently replicated spatial age effect is supplied by this reanalysis. Confidence intervals are exploratory within-section block intervals, not multiplicity-adjusted biological validation.
