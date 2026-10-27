# Muscle, matrix and vascular development in the postnatal heart

This is the consolidated rewrite of the biological synthesis and alternative fibroblast identity reports. It uses the clearer biological language discussed in chat, integrates the fibroblast findings into the developmental story, and keeps detailed analysis rules in Methods.

The full report has now been revised throughout using the preferred coronary–cardiomyocyte passage as its model. The abstract, introduction, results, discussion, figure legends, evidence tables and cell-boundary explanation describe the observation, its possible biological meaning and the remaining uncertainty. Methods retain the parameters needed to reproduce the work. “Lower/higher expression” describes changes in gene-expression scores without implying that fewer/more distinct genes were detected.

- [Rewritten combined Word report](Postnatal_heart_integrated_clear_language_report.docx)
- [Full report PDF](Postnatal_heart_integrated_clear_language_report.pdf)
- [All seven main and two supporting figure plates](Heart_atlas_integrated_figure_plates.pdf)
- [Readable text and figure legends](Report_text.md)
- [Checks of the numerical claims](data_claim_checks.json)
- [Figure inventory and original source paths](figure_inventory.json)

## What is incorporated

The report explains within-subtype CM maturation, later fibroblast matrix specialization, broader candidate fibroblast maps, the distinction between a narrow Gfpt2 identity and its more widely shared marker pattern, regional coronary expansion, earlier CM–non-CM relationships, and cross-age neighborhood comparisons. Inner myocardium is defined as muscle toward the chamber lining. The phrase “outside that envelope” is replaced by ordinary descriptions of distance from arterial structures.

The current alternative classifiers provide candidate labels for 7,992/10,867 eligible P21 fibroblasts. Direct scoring identifies the Gfpt2/Pi16/Adgrd1 pattern in 1,339 profiles, including 327 also showing the Col8a1 pattern. The report distinguishes these program-positive profiles from smaller, more consistently assigned subtype subsets, and includes the naming-marker-withheld evidence.

## Analysis status

This is an integrated writing and figure-presentation revision. No new biological model is fitted. New fibroblast classifier labels update identity and marker-pattern maps; the earlier CM association models, matched arterial-proximity comparison and transport plans have not been refitted using these labels. This distinction appears on the cover, in Results, Methods and the evidence table.

Original data and earlier releases are preserved:

- [Alternative fibroblast methods, annotations and source data](../fibroblast_multimethod_identity/README.md)
- [Earlier biological synthesis and source sensitivities](../cardiac_biological_synthesis/README.md)
- [Original detailed segmentation assessment](../cardiac_biological_synthesis/Heart_Atlas_segmentation_assessment.md)
- [Earlier subtype-conditioned CM models and transport](../cardiac_subtype_niches/README.md)

The matrix figure is rearranged from the existing validated tables. Other figures refer directly to the original PNG/PDF/SVG files recorded in the inventory, avoiding redundant copies. Both uploaded documents were inspected as source material; their text matched the prior releases, and the synthesis figure pixels also matched. Input-document fingerprints are recorded without redistributing the uploads.

## Reproduce the report

Run with a document Python providing python-docx, Pillow and PyMuPDF; the separate scientific Python provides matplotlib, pandas and NumPy for the rearranged matrix plate. LibreOffice renders the document for inspection.

```bash
python analysis/finalize_clear_heart_atlas_report.py --build
python analysis/finalize_clear_heart_atlas_report.py --resume
```

The source prose is [analysis/heart_atlas_clear_language_revision.md](../../analysis/heart_atlas_clear_language_revision.md). Validation checks numerical statements against the earlier source tables, current image hashes in Word, figure/PDF resolution, page boundaries, source links and file checksums. Earlier biological/model validation remains attached to the source analyses. One available spatial section per age and no new functional or lineage evidence remain material limitations.

## Attribution

Original atlas: Wang et al., Nature Cardiovascular Research (2026), https://doi.org/10.1038/s44161-026-00861-z . Vascular context: Lu et al., Nature Cell Biology (2021), https://doi.org/10.1038/s41556-021-00747-1 . Wormhole: Haviv et al., ICML (2024), https://arxiv.org/abs/2404.09411 . Public data derivatives retain attribution; journal PDFs are not redistributed.
