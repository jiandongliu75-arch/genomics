# Heart-atlas methods for relating CM development to other cardiac cells

- [Five-page Word assessment](Heart_atlas_nonCM_signaling_methods_review.docx)
- [PDF assessment](Heart_atlas_nonCM_signaling_methods_review.pdf)
- [Readable assessment](Report_text.md)
- [Measured signaling-target feasibility](signaling_target_feasibility.tsv)
- [Published WGCNA module coverage](published_WGCNA_module_coverage.tsv)
- [Proposed analysis specification](proposed_analysis_specification.json)
- [Validation](validation.json)

The paper combines spatial proximity, COMMOT ligand–receptor inference, random-forest/SHAP prediction, regulon inference and selected experimental follow-up. The useful next extension is source-resolved ligand–receptor–response modeling. Our previous identity/program Ridge models did not test the same features, learner or exact outcome implementation as the authors' model; their performance should not be treated as a direct replication comparison.

The review verifies the released 80:20 random cell split and recommends spatially separated validation, a CM receptor baseline, explicit CM-to-CM inputs, source-specific non-CM ligand features, and disjoint source neighborhoods. It documents assay limitations: Meis1 and Erbb2 are not measured targets, and two published RNA modules have only 2 and 10 measured members. Panel membership is not transcript detection, a verified ligand–receptor mechanism or pathway activation.

This is a **method review and feasibility audit**, not a newly fitted communication model or new evidence of a non-CM influence. It preserves all prior reports and observations.

Rebuild using the document Python environment:

```bash
python analysis/build_atlas_nonCM_signaling_review.py
soffice --headless --convert-to pdf --outdir results/atlas_nonCM_signaling_review \
  results/atlas_nonCM_signaling_review/Heart_atlas_nonCM_signaling_methods_review.docx
python analysis/build_atlas_nonCM_signaling_review.py --validate
```

The checked author-code commit is `102d80e6ea5de44f594461c727a4b3c6e1bc3508`. The paper Methods and released examples differ in database choice; exact original score and communication feature matrices are not recovered. See the assessment for evidence scope and biological interpretation.
