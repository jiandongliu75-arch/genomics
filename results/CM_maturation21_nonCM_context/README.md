# CM maturation and the developing cellular environment

Directly measured Xenium 21-gene CM scores, current non-CM maps, capillary-source sensitivity and fibroblast developmental interpretation. Original data and earlier reports are retained.

- [Biological Word report](CM_maturation_and_fibroblast_development.docx)
- [Full report PDF](CM_maturation_and_fibroblast_development.pdf)
- [Six integrated figure plates](CM_maturation_integrated_figures.pdf)
- [Source tables](CM_maturation_source_tables.zip)
- [Readable report with figures](Report_text.md)
- [Scoring definition](scoring_methods.json)
- [Analysis validation](analysis_validation.json) and [independent validation](independent_validation.json)

At P7, the full 21-gene comparison does not clearly establish higher overall maturation alongside capillaries. The lower-cycling observation from the original source definition is not reproduced with stricter capillary sources. Component-level relationships, transcript depth and source definition matter. These are observations from one spatial section per age, not causal or animal-level validation.

Col8a1-associated fibroblasts remain one family. Gfpt2/Pi16 is a descriptive name for the existing multi-gene candidate family; no double-positive gate or new classifier was introduced. Late Gfpt2 assignment disagreement remains unresolved. Expression similarities nominate developmental connections but do not establish fibroblast ancestry, migration or replacement. The report explicitly answers the questions about wall growth, early fibronectin, later collagen programs, Hortells, CM-derived TGFβ and the value of finer annotation.

## Per-cell data

Stage-qualified identifiers are required across sections. `cell_id` is unique within a section; `stage_cell_id` is globally stage-qualified. Coordinates are native Xenium micrometers.

| Stage | CM score table | CM/non-CM neighborhood table |
|---|---|---|
| P0 | [Scores](P0/CM_maturation21_scores.tsv.gz) | [Neighborhoods](P0/CM_nonCM_neighborhoods.tsv.gz) |
| P7 | [Scores](P7/CM_maturation21_scores.tsv.gz) | [Neighborhoods](P7/CM_nonCM_neighborhoods.tsv.gz) |
| P14 | [Scores](P14/CM_maturation21_scores.tsv.gz) | [Neighborhoods](P14/CM_nonCM_neighborhoods.tsv.gz) |
| P21 | [Scores](P21/CM_maturation21_scores.tsv.gz) | [Neighborhoods](P21/CM_nonCM_neighborhoods.tsv.gz) |

The score tables retain all previously scored CM candidates. Use `strict_ventricular_map_cohort` to recover the plotted ventricular selection. Models additionally use the eligibility conditions reported in Methods. `core_Z__maturation21` is the primary direct score; native, component, gene-group omission and 50-molecule sensitivity scores are explicitly named. Within-stage/identity percentiles are ranks, not fractions of physiological maturity.

Neighborhood fractions beginning `FB__` use all eligible fibroblast neighbors as denominator; `EC__` uses all eligible EC neighbors; `Broad__` uses all resolved neighbors. Unknown subtypes stay in lineage denominators. Consult [source coverage](nonCM_source_coverage.tsv) and [fibroblast assessability](fibroblast_assignment_coverage.tsv) before interpreting a zero: a family absent from the same-age reference is not assessable, not biologically absent. Centroid distance columns are separate descriptive measurements; they do not establish contact.

## Main evidence tables

- [Matched capillary results](matched_capillary_maturation_results.tsv): original pairs, strict ventricular pairs, identical depth-eligible pairs, and newly fitted QC-source pairs. Twenty strict ventricular core comparisons form the separately corrected matched family. Other rows are sensitivity analyses.
- [QC-source matched pairs](QC_capillary_refitted_pairs.tsv.gz): newly fitted without using outcome scores.
- [Non-CM associations](nonCM_maturation_associations.tsv): 235 supported primary tests, corrected together. Additional CM-state, RNA-view, radius and score variants are sensitivities. Primary points are conditional associations, not direct cell influences.
- [Incremental prediction](incremental_spatial_prediction.tsv): capture/size/CM-identity baseline, anatomy/broad context, added fine identities, and added measured fibroblast programs. Five buffered spatial folds are not independent hearts.
- [Single-gene omissions](leave_one_gene_out_correlations.tsv), [equal-depth checks](depth_matching_checks.tsv), [same-cell stage means](same_CM_depth_and_gene_omission_stage_means.tsv), and [additive gene contributions](gene_contributions_to_stage_score_change.tsv).
- [Fibroblast matrix changes](fibroblast_matrix_gene_stage_changes.tsv), [library matrix scores](fibroblast_library_matrix_scores.tsv), [expression continuity candidates](fibroblast_candidate_expression_continuities.tsv), and [spatial geometry snapshots](fibroblast_spatial_geometry_snapshots.tsv).

The score is a documented frozen-reference standardized mean using the final paper's 21 measured genes. It is not asserted to reproduce the authors' unavailable precomputed Seurat AddModuleScore object. Controls and reference scaling are not reselected during omission tests. Fixed-depth sampling removes information and can remove biological capture/size differences; it does not estimate a fraction of technical artifact.

## Reproduction

From the repository root, with the existing scientific Python environment and frozen input files:

```bash
OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 \
  /workspace/heart-atlas/venv/bin/python analysis/assess_CM_maturation21_neighborhoods.py
OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 \
  /workspace/heart-atlas/venv/bin/python analysis/validate_CM_maturation21_context.py
OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 \
  /workspace/heart-atlas/venv/bin/python analysis/plot_CM_maturation21_context.py
python analysis/write_CM_maturation21_report.py
soffice --headless --convert-to pdf --outdir results/CM_maturation21_nonCM_context \
  results/CM_maturation21_nonCM_context/CM_maturation_and_fibroblast_development.docx
python analysis/finalize_CM_maturation21_report.py
```

After a completed release, `python analysis/finalize_CM_maturation21_report.py --resume` verifies source and artifact hashes without refitting. Figure labels, score reconstruction, independent neighbor counts, regression coefficients, QC matching constraints, corrected testing families, buffered folds, Word image hashes and rendered layout are checked. No new fibroblast fate trajectory, segmentation, transport model or biological replicate is claimed.
