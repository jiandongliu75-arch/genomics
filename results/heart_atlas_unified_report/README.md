# Muscle, matrix and vascular development in the postnatal heart

This is the single integrated report containing the CM, endothelial, fibroblast, vascular and cross-age neighborhood story, with the latest completed fibroblast–CM refit incorporated. The narrative is rewritten as one account; superseded conclusions are replaced rather than appended as a separate report.

- [Complete Word report](Postnatal_heart_unified_biological_report.docx)
- [Complete PDF report](Postnatal_heart_unified_biological_report.pdf)
- [All nine main and four supporting figure plates](Postnatal_heart_unified_figure_plates.pdf)
- [Consolidated source-table download](Postnatal_heart_unified_source_tables.zip)
- [Readable narrative and figure legends](Report_text.md)
- [Source-table index](source_table_index.tsv)
- [Figure inventory and original source paths](figure_inventory.json)
- [Document and release checks](release_validation.json)

## What is integrated

The report covers broad cardiac populations and the five practical EC categories; ventricular CM identity and within-subtype maturation; full-RNA fibroblast matrix features and named P21 substates; current candidate maps, shared patterns and Gfpt2 uncertainty; published perinatal coronary growth; earlier capillary/CM comparisons; the latest FB–CM/endocardial refit and boundary sensitivities; calibrated exact/Wormhole correspondence with updated fibroblast descriptions; and a focused experimental roadmap.

The latest refit supplies the current conclusions. The earlier P7 Postn metabolic association weakens with broader assignments. No FB identity/program association passes the full 925-test primary correction. A P0 elastic-matrix/sarcomere relationship remains a follow-up candidate. The P14 Ankrd1/endocardial structural association is the clearest corrected relationship but depends on transcript assignment. Broader P14 Col8a1 candidates are more often near arterial structures than matched Gfpt2 candidates, with a smaller contrast, approximately 11% versus 7%.

Absent same-age reference classes are marked **not assessable**, rather than biological zeros. Subtype identity, a shared gene pattern and developmental fate remain separate. The measured 21-gene CM panel and all seven FB program definitions are included. Kcnd2/Irx5 were not measured in the deposited 403-gene Xenium panel.

## Analysis status

This release assembles existing validated results and four new presentation plots; no new biological model is fitted. The latest refit is incorporated from [the completed FB–CM release](../cardiac_fibroblast_CM_refit/README.md). Earlier capillary-adjacency and nuclear cycling comparisons retain their original cohorts and uncertainty definitions. Original exact/Wormhole plans remain frozen; their fibroblast descriptions are updated. Source tables identify which version each result uses.

Original individual reports, processed matrices and source annotations are preserved. Figure files already available elsewhere are referenced directly instead of redundantly copied. The new matrix plate replaces the earlier stringent arterial-comparison panel with the latest broader-candidate contrast.

## Source data

The consolidated ZIP contains the latest FB–CM source tables under `fibroblast_CM_refit/`, plus relevant full-RNA matrix/marker tables, alternate-identity audit tables, earlier vascular contrasts and exact/Wormhole calibration tables. It does not duplicate raw counts or the large per-CM NPZ arrays. These arrays remain in the completed refit folder and their schema is described there. The source-table index records original paths and hashes.

## Reproduce the document

Use document Python with python-docx, Pillow and PyMuPDF; scientific Python supplies NumPy, pandas, SciPy and matplotlib. LibreOffice renders the Word file.

```bash
python analysis/finalize_unified_heart_atlas_report.py --build
python analysis/finalize_unified_heart_atlas_report.py --resume
```

The editable prose is [analysis/heart_atlas_unified_report.md](../../analysis/heart_atlas_unified_report.md). Validation checks current figure hashes in Word, numerical claims, table/reference consistency, obsolete-status language, rendered pages, resolution, figure bounds, source links and checksums. Original biological-model validation remains in the source releases. One available section per age and the absence of direct functional/lineage evidence remain explicit limitations.

Original atlas: Wang et al., Nature Cardiovascular Research (2026), https://doi.org/10.1038/s44161-026-00861-z . Vascular context: Lu et al., Nature Cell Biology (2021), https://doi.org/10.1038/s41556-021-00747-1 . Wormhole: Haviv et al., ICML (2024), https://arxiv.org/abs/2404.09411 . Journal PDFs and user uploads are not redistributed.
