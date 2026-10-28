# Capillary proximity and CM development

[Biological report](Report_text.md) · [Word report](Capillary_maturation_retest.docx) · [PDF report](Capillary_maturation_retest.pdf) · [Integrated figure](Capillary_maturation_integrated_figure.pdf)

The new test supports a small adjusted proximity association at P7. It does not establish endothelial PTPRM or oxygen as the mechanism. Exact author score/masks/graph are unavailable; one spatial section per age limits replication. See the report and methods.json for primary and sensitivity definitions.

Reproduce with `OPENBLAS_NUM_THREADS=2 /workspace/heart-atlas/venv/bin/python analysis/retest_capillary_maturation.py`, then `analysis/report_capillary_maturation_retest.py` in that environment, followed by `python analysis/write_capillary_maturation_retest.py`. Existing counts, cell identities and boundaries are unchanged.
