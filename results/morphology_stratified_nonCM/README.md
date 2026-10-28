# Morphology, non-CM subtypes and CM development

[Word report](Morphology_nonCM_CM_development.docx) · [PDF report](Morphology_nonCM_CM_development.pdf) · [Narrative](Report_text.md) · [Morphology-only review viewer](Morphology_only_ROI_review.html) · [Frozen editable regions](morphology_ROIs.geojson)

The annotations are outcome-independent DAPI-based candidates. Definitive LV/RV, papillary and large-vessel identities remain unresolved. The usable compact/cavity-facing patches cover about 31% of strict CMs; no original author masks were recovered. The new subtype screen supports a P0 Postn-dominated fibroblast proximity association, without establishing subtype-specific causation. P7 capillary nearness remains positive but does not pass correction across the broader 48-test family.

All supported sources, alternative definitions, uncertain identities and unassessable tests are retained. Fine immune/mural/epicardial subtypes are not available and are not invented. The primary CM score is the frozen measured 21-gene score, not physiological maturity or exact author AddModuleScore.

Reproduce: annotate_morphology_regions.py; morphology_nonCM_associations.py; check_morphology_score_robustness.py; report_morphology_nonCM.py with the scientific Python environment; morphology_review_viewer.py and write_morphology_nonCM_report.py with base Python; validate_morphology_nonCM.py with scientific Python. acquire_anatomy_images.py is an optional retrieval helper: no antibody-channel download completed in this turn.
