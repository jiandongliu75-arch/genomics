# Combined snRNA-seq and Xenium analysis of CM development

This analysis tests whether source-specific molecular inputs from surrounding cardiac cells improve prediction of directly measured CM developmental expression after receiver state, anatomy, RNA capture and CM-to-CM inputs are considered.

- [Editable biological report (Word)](Joint_RNA_Xenium_CM_development.docx)
- [Biological report (PDF)](Joint_RNA_Xenium_CM_development.pdf)
- [Four integrated figure plates (PDF)](Joint_RNA_Xenium_integrated_figures.pdf)
- [Report text and individual figures](Report_text.md)
- [Analysis cohorts](analysis_cohorts.tsv)
- [Nested-model prediction gains, intervals and primary q values](incremental_prediction_gains.tsv)
- [All held-out performance results, including CM states and both learners](OOF_performance.tsv)
- [Source-feature sensitivities](source_feature_sensitivity_summary.tsv)
- [Biologically selected RNA availability](biologically_selected_RNA_availability.tsv)
- [Ankrd1 versus working-CM independent RNA programs](RNA_Ankrd1_working_independent_program_contrasts.tsv)
- [Full RNA source–pair support](RNA_pair_source_support.tsv.gz)
- [Selected candidate pairs](selected_RNA_supported_pairs.tsv)
- [Non-conduction CM receptor sensitivity](RNA_selected_pair_receiver_sensitivity.tsv)
- [Source-table index and checksums](source_table_index.tsv)
- [Independent RNA validation](RNA_independent_validation.json), [independent spatial validation](spatial_independent_validation.json) and [release validation](release_validation.json)

The biological conclusions and their limits are given in the report. A positive prediction gain does not indicate a positive biological effect on maturation. A negative gain is worse prediction, not inhibition by a cell source. RNA-supported compatibility is not activated signaling, perfusion or deposited matrix.

## What is measured and inferred

The source reference contains 73,052 previously QC-filtered, DecontX-corrected nuclei, with two RNA libraries at each of P0, P7, P14 and P21. All 180 selected age–pair qualifications retain receptor support after excluding provisional conduction and uncharacterized CM candidates. These libraries do not create additional spatial-heart replicates.

The four Xenium sections contain 403 measured targets. The CM 21-gene and component outcomes remain directly measured. No unmeasured RNA program is inserted into Xenium. Core denotes the existing high-confidence transcript-assignment view, rather than strictly nuclear counts. Core and deposited-RNA models use the same cells. CMs lacking a measured nucleus are not automatically excluded.

The primary database is mouse CellChat packaged with COMMOT 0.0.3. CellPhoneDB v4 pair agreement is annotated. Candidate eligibility requires all subunits, replicated RNA availability and no primary outcome-gene overlap. Fine arterial/lymphatic source-specific molecular tests lack adequate RNA-library support; observed source counts and pooled endothelial inputs still retain those profiles. Unknown or unavailable families are not biological zeros.

## Cell-level data

Each age directory (`P0`, `P7`, `P14`, `P21`) contains:

- `core_CM_candidate_map_values.tsv.gz` and `native_CM_candidate_map_values.tsv.gz`: CM IDs, coordinates, identity, domain, test eligibility, five observed CM outcomes, selected molecular candidates and held-out forest predictions. NA means unavailable or not predicted; it is not zero. Candidate input columns are log1p model-allocated inputs, not receptor activation.
- `core_OOF_predictions.npz` and `native_OOF_predictions.npz`: all held-out predictions, outcomes, eligibility, identities, learner names and feature-set names. Predictions outside the eligible cohort are NaN.
- `source_domains.tsv.gz` and `receiver_domains.tsv.gz`: recorded cells, positions, source identities and domain eligibility.
- `feature_definitions.tsv` and `stage_pairs.tsv`: molecular feature definitions, source groups, ranges and stage-specific pairs.
- Per-view numerical transport diagnostics, source–pair aggregation, fold performance, fold checks and source-feature perturbations.

Larger `*_predictor_features.npz` matrices are regenerable local intermediates and are not included in Git. Full per-edge transport plans are aggregated by source identity before export and are not distributed as a separate edge atlas. Original count matrices and previous QC outputs remain in the upstream repository releases and recorded workspace inputs.

## Reproduction

Use the scientific Python environment with NumPy, pandas, SciPy, h5py, anndata, scikit-learn and Matplotlib. Word/PDF writing uses python-docx, Pillow, ReportLab and PyMuPDF in the document environment. Package versions are recorded in `analysis_environment.json`.

```bash
python analysis/joint_signaling_RNA.py
python analysis/joint_signaling_spatial.py --solver-check
python analysis/joint_signaling_spatial.py --resume
python analysis/joint_signaling_models.py --resume
python analysis/joint_signaling_models.py --summarize
python analysis/validate_joint_signaling.py --RNA-only
python analysis/validate_joint_signaling.py
python analysis/export_joint_signaling.py
python analysis/plot_joint_signaling.py
python analysis/compose_joint_signaling_report.py
python analysis/write_joint_signaling_report.py
python analysis/finalize_joint_signaling.py
```

Scripts resolve repository inputs relative to the checkout. The full corrected RNA H5AD is expected at `/workspace/heart-atlas/results/cell_identities_subtypes/snrna_qc_annotated_subclustered.h5ad`, as documented by the existing QC pipeline; adapt that path if relocating the analysis. Raw source accessions are GSE290989 and Zenodo records 14961929/14961942. [Input/code checksums](source_checksums.tsv), [released-file checksums](file_manifest.tsv) and [configuration](configuration.json) identify this analysis. Numerical optimizers and parallel forests can have small platform-dependent differences within the recorded convergence tolerances.

Four spatially separated held-out domains per section have a 200-µm interior receiver margin and at least a 400-µm train/test receiver gap. All source types compete jointly inside each independent domain. This uses COMMOT's collective COT component, not the default four-solver ensemble or the authors' recovered feature matrix. Native/core, learner and CM-state results are supplied as sensitivities. Conditional spatial-block intervals do not establish independent-heart replication or causality.

The linear sensitivity uses a minimum scale of 0.05 for predictors beyond baseline A to prevent numerical extrapolation from near-zero training variation. This adjustment does not change the primary random-forest predictions. The report records all primary comparisons rather than selecting a favorable pathway or age.
