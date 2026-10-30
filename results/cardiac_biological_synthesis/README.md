# Regional coordination of muscle, matrix and vascular development

This revision places the biological questions first: how ventricular CM identities and maturation states relate to fibroblast matrix specialization, regional vascular expansion and later vascular integration. It synthesizes Lu et al. (2021) without treating its signaling pathway as the main explanation for our data.

## Read the revised report

- [Editable Word manuscript](Postnatal_heart_biological_synthesis.docx)
- [Rendered manuscript PDF](Postnatal_heart_biological_synthesis.pdf)
- [Six integrated main figures](Biological_synthesis_integrated_figures.pdf)
- [Manuscript text](Manuscript_text.md)
- [Dedicated assessment of segmentation in the Heart Atlas](Heart_Atlas_segmentation_assessment.md)
- [Figure source tables and candidate annotations](Biological_synthesis_source_tables.zip)
- [Publication figures: PNG, vector PDF and editable SVG](Biological_synthesis_publication_figures.zip)

## Biological conclusions and their scope

- CM identity and maturation vary independently within early sections.
- Early Postn-prominent fibroblast remodeling is followed by a more specialized RNA landscape. Residual late Postn detection is distinguished from persistence of a distinct Postn family.
- Col8a1 fibroblasts have a stronger late basement-matrix-associated RNA program than Gfpt2 fibroblasts, with a partial vascular-associated component in the earlier confidence-supported comparison.
- Published perinatal vascular work establishes regionally concentrated expansion of the existing coronary plexus, especially into inner myocardium; endocardium is distinct from that source.
- EC/CM cycling, capillary apposition and later CM metabolic/structural specialization are different measurements. We do not equate a nearby EC with measured perfusion or a single maturation mechanism.
- The P7 cycling-CM/Postn-FB metabolic association weakens with broader FB evidence. The P14 Ankrd1-CM/endocardial-EC sarcomere association persists in this sensitivity. Both remain observational and geometry-sensitive.
- Transport and RNA resemblance nominate developmental comparisons, not cell ancestry or transition rates.

## Fibroblast reassessment

Frozen atlas labels are not overwritten. Nuclear markers are optional corroboration in new candidate layers; supported alternative nuclear calls remain counterevidence. Whole-cell and boundary-core calls must agree. The additional multi-gene layer removes a compulsory individual naming-marker anchor but retains at least two supplied markers, reference-library/raw-reference consistency, weight/margin thresholds and capture/source QC in both cell and core views. These operational candidates are not independently validated identities or calibrated posterior probabilities.

| Stage | Legacy strict sources | Nuclear-optional + anchor | Multi-gene boundary candidates | Eligible broad FBs |
|---|---:|---:|---:|---:|
| P0 | 1,120 | 1,240 | 1,249 | 2,020 |
| P7 | 3,656 | 4,179 | 4,185 | 7,365 |
| P14 | 589 | 795 | 2,884 | 8,098 |
| P21 | 679 | 941 | 3,707 | 10,867 |

Removing nuclear corroboration alone does not resolve late Gfpt2 placement. Among 2,716 eligible public P21 Gfpt2-labelled FBs, only 129 receive a Gfpt2 whole-cell proposal before nuclear filtering. The broader evidence layer retains that disagreement rather than forcing a label or concluding biological absence. Full [coverage](FB_assignment_coverage_comparison.tsv), [gate audit](FB_assignment_gate_audit.tsv), [proposal cross-tabulation](FB_public_vs_whole_cell_proposals.tsv) and per-cell evidence are provided.

The primary analysis remains available in [the preceding release](../cardiac_subtype_niches/README.md). New component models are [exported in full](nuclear_optional_FB_maturation_coefficients.tsv), including 30-µm/contact sensitivities. Their denominators retain all eligible broad FBs, including unresolved families. Revised maps and earlier maps are explicitly versioned; the transport representation/plans are unchanged.

## Segmentation assessment

The authors used membrane/interior/nuclear staining and manually benchmarked 10 fields against nucleus expansion. Supplementary Figure 1 supports improvement, but shows modest absolute boundary agreement and marked deterioration at P21 at the demanding 90% overlap threshold. Approximate visual scores are not recovered raw benchmark values. The object score does not measure the fraction of transcripts correctly assigned. The full [assessment](Heart_Atlas_segmentation_assessment.md) separates what the paper demonstrates from what remains uncertain at CM–non-CM interfaces. The downloaded journal supplement is not redistributed.

## Reproduce

Use the existing scientific environment for analysis/figures and a Python with python-docx, Pillow and PyMuPDF for Word/release validation. LibreOffice renders Word for layout checks.

```bash
OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 \
  /workspace/heart-atlas/venv/bin/python analysis/run_cardiac_biological_synthesis.py --document-python python
```

The source prose is [analysis/cardiac_biological_story.md](../../analysis/cardiac_biological_story.md). The segmentation assessment is an editorial review linked to the cited article Methods and directly inspected Supplementary Figure 1. Scientific validation includes independently reconstructed neighbor fractions and centered least-squares refits; document validation checks current figure hashes, resolution, page boundaries and source links. One available section per age and no new biological/functional validation remain material limitations.

## Attribution

Original atlas: Wang et al., Nature Cardiovascular Research (2026), https://doi.org/10.1038/s44161-026-00861-z . Vascular context: Lu et al., Nature Cell Biology (2021), https://doi.org/10.1038/s41556-021-00747-1 . Wormhole: Haviv et al., ICML (2024), https://arxiv.org/abs/2404.09411 . Public data derivatives retain attribution; uploaded journal PDFs are not included.
