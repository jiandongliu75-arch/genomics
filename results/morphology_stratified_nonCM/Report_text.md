# Cardiac anatomy, neighboring cell subtypes and CM transcriptional development

## What the new analysis adds

The new analysis separates recognizable wall and cavity-facing regions before asking how nearby non-CMs relate to CM development. It retains a small positive P7 capillary estimate and identifies a P0 fibroblast association in conservative compact-wall patches. The relationships change with age, CM program and neighborhood definition. They do not support a single rule that one neighboring cell type uniformly accelerates CM maturation.

Anatomical certainty is a material limit. I annotated the available DAPI images without displaying maturation scores, RNA expression or cell-type maps during ROI construction, and froze the resulting polygons before the new model fits. The annotations identify conservative continuous compact-wall patches and visible cavity-facing myocardium. They do not establish definitive left/right ventricular masks or verified papillary muscles. Cavities are partly collapsed or obliquely cut, and the acquired antibody channels are no longer cached. This is an outcome-independent candidate annotation, not an independently blinded expert histology assessment or recovery of the authors' masks.

## What can be recognized from morphology?

There are 31 candidate ROIs across the four sections: 16 compact-wall patches, seven cavity-facing patches, two possible papillary/trabecular projections and six vessel-like surroundings. The latter two categories have low confidence and are excluded from primary tests. Only the compact and cavity-facing patches enter the main anatomical analysis (Figure 1). They contain 36,647 frozen ventricular CM profiles, about 31% of the original strict cohort; the remaining profiles are retained as unclassified or low-confidence rather than forced into anatomical categories.

| Age | Original strict CMs | Compact-wall CMs | Cavity-facing CMs |
|---|---:|---:|---:|
| P0 | 17,272 | 5,529 | 0 |
| P7 | 30,485 | 10,372 | 1,600 |
| P14 | 33,659 | 8,890 | 1,693 |
| P21 | 37,949 | 7,776 | 787 |

The cavity-facing label means visible chamber-facing or ridged myocardium. It should not be treated as a verified histological trabecular layer. At P0, no reliable cavity-facing patch was designated. Image-coordinate territories are retained for review, but anatomical LV/RV names remain unresolved. The possible papillary projections cannot be distinguished confidently from trabecular ridges in these two-dimensional images. Closed, DAPI-poor lumens define vessel-like candidates, but do not establish artery or vein identity. These small candidate regions also contain too few independent spatial clusters for reliable region-specific association tests.

A self-contained morphology-only viewer accompanies the report. It displays the DAPI images and ROI outlines without RNA or maturation scores, supports vertex editing and saves reviewed polygons as GeoJSON. This makes the anatomical interpretation concrete and reviewable. The current statistical results use the frozen version, not any subsequently edited copy.

## Fibroblasts: the strongest early relationship is broader than a subtype-specific effect

At P0, CMs closer to fibroblasts in the annotated compact-wall patches have higher 21-gene scores. The Postn-associated estimate is +0.048 (+0.026 to +0.070) SD, with q=0.0043 across the 48 primary comparisons. The metabolic-regulatory component is also higher. Deposited RNA and cell-center distances give the same direction. Targeted score checks retain the positive association after omitting Myl2, jointly omitting Cebpb/Rxrg/Lpin1, and equalizing measurement to 50 molecules per CM. This reduces concern that a single panel gene or capture depth entirely creates the result.

However, 1,946 of the 2,020 eligible P0 fibroblasts are Postn-associated. The broad-fibroblast and Postn estimates are therefore nearly identical. This is evidence for a relationship with the early Postn-dominated fibroblast environment; it cannot establish that Postn fibroblasts have a unique effect compared with other P0 fibroblast subtypes. It also does not demonstrate that fibroblasts cause CM maturation. Local matrix organization and shared regional development remain possible explanations.

At P7, the overall-score relationships with Postn-associated and Gfpt2/Pi16-associated fibroblast proximity are uncertain. Later Col8a1-associated and Gfpt2/Pi16-associated proximity also fails to provide a corrected overall-score association. Lsamp-associated candidates are tested where spatial support permits; Kirrel3 candidates lack adequate local support in these annotated patches. No general maturation-promoting role is established. Col8a1 remains one family, and Gfpt2/Pi16 is a multi-gene interpretation rather than a mandatory double-positive gate. Missing same-age reference families and inadequate local coverage are recorded as not assessable, not biological absence.

Measured fibroblast programs provide a complementary result. At P7, CMs surrounded by stronger Eln/Fbln1/Fbln2 expression have slightly lower overall maturation scores: -0.028 (-0.044 to -0.012), secondary q=0.047. The direction persists in deposited CM RNA and after controlling for informative fibroblast count and RNA capture. This is an exploratory matrix-context relationship. The genes do not measure elastin deposition, tissue stiffness or an outflow-tract location, and the finding does not identify a particular fibroblast family as causal.

## Endothelial cells: vascular proximity remains plausible, but maturation is not one uniform response

Within the annotated P7 regions, the capillary estimate remains positive: +0.020 (+0.004 to +0.036) SD. Its size is similar to the preceding within-region retest. It is less precise in this smaller anatomical subset, and q=0.119 after correcting across the full 48-comparison source screen. Deposited RNA has a positive estimate with an interval crossing zero. A larger fraction of local ECs being capillary-like does not show the same clear relationship. Nearness and relative neighborhood composition ask different questions; neither should be substituted for the other to claim a stronger result.

Thus, adding morphology has not reversed the P7 capillary estimate, but the new screen does not independently confirm it across multiple animals or establish a universal vascular effect. Gene-omission checks retain a positive estimate, while equal-molecule scoring in this smaller anatomical subset gives a positive but uncertain estimate. The earlier targeted analysis and this broader screen also use different cohorts and multiplicity families. A change in statistical significance alone is not evidence that regional anatomy explains the association.

P14 venous-like EC proximity gives another candidate: +0.029 (+0.009 to +0.050), q=0.075. It does not pass primary correction, and the local venous fraction is not clearly associated with the overall score. Arterial-like, endocardial-like and cycling-associated EC proximity do not provide corrected overall-score findings in the primary screen. Lymphatic-like cells and some stage-specific EC classes lack sufficient local support and cannot be assessed reliably.

Endocardial-like ECs illustrate why the individual CM programs matter. In the annotated P14 cohort, nearness to these cells relates to lower sarcomere RNA, including within Ankrd1-associated CMs. This differs from the earlier positive association using an endocardial neighborhood-composition measure in a broader cohort. The present fraction-based comparison is uncertain. The earlier positive relationship should therefore not be treated as an established general effect of endocardial ECs on CM structural maturation. Distance, composition, regional selection and technical definitions need to be reconciled before assigning a biological direction to that relationship. Neither result establishes greater overall CM maturity.

## Mural, immune and epicardial populations

Pericytes, smooth-muscle cells, immune cells and epicardial cells are included as supported broad sources. None provides an overall-score association surviving the primary correction. Some component-specific or region-specific estimates differ, but they remain secondary hypotheses in one section per age. The current spatial reference does not resolve validated fine immune, pericyte, smooth-muscle or epicardial subtypes, so those labels have not been invented or inferred from location alone.

The most useful biological picture is therefore a changing set of local environments. Early fibroblast-rich matrix settings, capillary proximity and particular matrix-expression programs accompany different parts of CM transcriptional development. The data identify comparisons to test, rather than a ranked list of cells that promote or inhibit maturation. Independent hearts, expert-reviewed regional anatomy and functional CM measurements are needed to establish which associations extend beyond these sections.

## Methods and interpretation limits

The ROI construction script reads only cached DAPI overviews and physical image extents. Frozen GeoJSON and image checksums record the exact annotations. Boundary and alphaSMA/vimentin TIFF caches were absent; Zenodo retrieval failed with a network proxy403 response even after network access was granted. Consequently, no newly recovered antibody images or verified vessel walls are claimed. Previously validated RNA counts, segmentation and cell identities were retained, including valid anucleate CM sections.

The CM outcome is the previously frozen mean of 21 measured Xenium gene Z values, with metabolic-regulatory, sarcomere, excitation–contraction and cycling components separate. It is not the authors' unavailable exact AddModuleScore object, and it does not measure physiological function. Separate deposited and core RNA outcomes use fixed full-stage scales. The primary source metric is minus log1p shortest weighted centroid-path distance on the frozen tissue-screened Delaunay graph; paths are not blood flow, diffusion or verified membrane contact. Cell-center nearness and 50-µm source fractions are retained as alternative questions.

Models compare CMs within image ROI, CM subtype and 200-µm grid, with flexible capture and cell-size controls, missing nuclear area and continuous cavity/surface distances. A source requires at least 30 eligible profiles, at least 100 locally supported CMs, and support in at least eight 400-µm clusters; each fitted stratum also requires adequate local source support across eight clusters. Insufficiently supported subtype/region tests are explicitly skipped. Intervals use clustered CR1 t uncertainty within the observed section. They do not incorporate independent-heart, annotation or segmentation variability.

The 48 primary tests cover supported broad sources, current EC subtypes and current fibroblast families across the four ages, using overall core CM score in the annotated regions. Broad parents and subtypes can overlap substantially; a positive subtype result is not proof of a unique effect beyond its parent lineage. All core component, alternative metric, source-rule, CM-state and program comparisons are corrected together as secondary tests; native RNA is retained as a measurement sensitivity. The measured fibroblast-program model includes informative-source count and capture; a version without those source controls is also retained. Nuclear-supported and public-plus-broad-supported EC source rules are distinguished, and unresolved identities remain in lineage denominators.

The named P0 Postn and P7 capillary comparisons receive targeted follow-up score checks, including Myl2 omission, joint Cebpb/Rxrg/Lpin1 omission and equal-50-molecule scores in the same depth-eligible cells. These 10 sensitivity tests are corrected separately and do not add primary discoveries or independent biological replication.

No new biological samples, segmentation, spatial gene imputation, lineage tracing or causal perturbation were added. Chamber, papillary and vessel-like candidates require anatomical review before they can support definitive compartment-specific conclusions. The exact original author masks remain unrecovered.

Source: Wang et al., Spatially guided in vivo single-cell functional genomics of postnatal heart, Nature Cardiovascular Research, DOI 10.1038/s44161-026-00861-z. https://pmc.ncbi.nlm.nih.gov/articles/PMC13561863/ The deposited source datasets are GSE290989 and Zenodo 14961929/14961942; the reviewed reference counts, source definitions and earlier QC remain unchanged.
