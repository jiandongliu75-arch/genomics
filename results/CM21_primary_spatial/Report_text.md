# The 21-gene CM maturation landscape and neighboring cardiac cells

## The question answered here

Where are ventricular cardiomyocytes with higher or lower 21-gene maturation scores at each postnatal stage, and which non-CM subtypes surround them? Every maturation map, comparison and model in this report uses the same 21 genes measured directly by Xenium. The existing RNA-supported annotations identify the neighboring cell populations.

The maps retain all 131,605 reviewed ventricular CM profiles passing the existing QC thresholds. Anatomical classification and agreement with public labels are not requirements for inclusion. CMs in uncertain anatomical regions and valid CM sections without a visible nucleus remain in the maps. Higher and lower groups are defined separately within each age, so a higher-scoring P0 CM is not assumed to have the same absolute maturity as a higher-scoring P21 CM.

| Stage | All QC-eligible ventricular CMs | Lower 21-score quartile | Higher 21-score quartile |
|---|---:|---:|---:|
| P0 | 17,981 | 4,495 | 4,496 |
| P7 | 35,112 | 8,778 | 8,779 |
| P14 | 37,740 | 9,435 | 9,436 |
| P21 | 40,772 | 10,193 | 10,194 |

## Where are higher- and lower-scoring CMs?

Figure 1 shows the directly measured 21-gene score on each section, using the same color scale across ages. Figure 2 highlights the higher and lower quartiles within each stage. Together, these maps show the age-related shift in the index and the remaining variation among CMs in the same section.

The complete cell-level tables contain each CM's coordinates,21-gene score, within-stage rank, high/low group and distances to supported non-CM sources. The interactive viewer allows selection of a stage and neighboring cell subtype, followed by inspection of individual CM scores and local source positions. All 21 targets are present in the Xenium panel; individual cells need not detect every target. The median number detected per CM is 10 at P0, 11 at P7 and 13 at P14/P21.

## Coronary capillary proximity: the raw and adjusted answers

At P14, higher-scoring CMs are closer to capillary ECs in the direct high-versus-low comparison, and the continuous-score analysis gives the same direction. At P7, the raw high-versus-low comparison does not show that pattern: the highest-scoring quartile is slightly farther from the QC-supported capillary sources. After accounting for RNA capture, cell size, CM identity and local position, a small positive capillary association appears. The statement that higher-scoring CMs are closer to capillaries at P7 is therefore conditional on those adjustments; it is not a simple unadjusted rule.

| Stage | Lower-score CMs: mean nearest capillary, µm | Higher-score CMs: mean nearest capillary, µm | Adjusted association, SD | All-source BH q |
|---|---:|---:|---:|---:|
| P0 | 53.6 | 51.6 | +0.019 | 0.165 |
| P7 | 44.2 | 46.4 | +0.028 | 5.83e-09 |
| P14 | 41.6 | 38.8 | +0.024 | 9.82e-06 |
| P21 | 53.9 | 51.5 | +0.008 | 0.0822 |

Positive adjusted coefficients mean that a higher 21-gene score accompanies greater capillary nearness. The adjusted coefficients are +0.028 at P7 and +0.024 at P14. These directions persist with deposited RNA and when anatomical labels are added while retaining unclassified regions. The P0 and P21 capillary comparisons do not pass the same all-source correction. This does not establish a biological onset or disappearance of a vascular relationship.

RNA capture is strongly related to the measured index and also varies spatially. That is why the raw and adjusted views are displayed together. Adjustment can reduce technical confounding, but it does not make one observational model a definitive account of maturation or establish that neighboring vessels cause the score differences.

## Which fibroblast subtypes accompany different 21-gene scores?

The full-cohort continuous analysis identifies a small positive relationship with Gfpt2/Pi16-associated fibroblast nearness at P7: coefficient +0.015, q=0.047. In contrast, Lsamp-associated fibroblast nearness relates to lower scores at P14 (-0.027, q=0.019) and P21 (-0.017, q=0.001). These directions persist in deposited CM RNA and the anatomy-covariate sensitivity.

These are conditional spatial associations. The raw high/low summaries do not always show the same contrast, so the results should not be translated into claims that these fibroblasts promote or suppress maturation. Nor does an Lsamp label establish an outflow-tract location.

Postn-associated fibroblast proximity does not show a corrected whole-cohort relationship at P0 or P7. Col8a1-associated fibroblast proximity does not show a corrected whole-cohort relationship at P14 or P21. Later Gfpt2/Pi16 assignments remain uncertain, and no corrected P21 association is established. A relationship in a selected compact-wall patch should not be generalized to the entire section. Unavailable reference families and inadequately supported sources remain not assessable, rather than being counted as biologically absent.

## Other neighboring cardiac cells

The five supported EC classes are evaluated where sufficient cells are available. Endocardial-like EC nearness relates to lower 21-gene scores at P7 in the adjusted model. Arterial-like and venous-like sources do not show corrected associations in this full-cohort screen. Lymphatic-like sources lack sufficient support for a reliable comparison. Broad EC results overlap with their capillary subtype results and are not independent evidence for a separate mechanism.

Among other broad populations, pericyte nearness relates to slightly higher scores at P14 and P21; immune-cell nearness does so at P14. Smooth-muscle nearness relates to lower scores at P21, and epicardial nearness relates to lower scores at P7–P21 in the adjusted models. These patterns can reflect tissue organization and regional CM states. They do not establish trophic, inhibitory or signaling effects, and validated fine immune, mural and epicardial subtypes are not supplied by the present reference.

Figure 3 presents every supported source rather than only selected findings. The left panel shows the raw score–nearness relationship, and the right panel shows the adjusted continuous-score relationship. The full high/low neighborhood table also gives distances and the proportion of CMs with each source within 50µm. These are alternative descriptions of the same 21-gene outcome, not additional maturation programs.

## What does strict high/low matching add?

The primary continuous analysis uses all eligible CMs. As a separate sensitivity, higher- and lower-score CMs were matched within CM identity and local regions for RNA capture, cell area and missing nuclear area, without matching on neighboring cells. Only 1,004 disjoint pairs across the four ages satisfy the strict calipers. None of the matched distance comparisons survives correction across the source screen.

This is important uncertainty. The matched analysis does not independently confirm the source-specific findings, and its small selected subset cannot replace the all-CM question. The strongest defensible interpretation is that the 21-gene maps show within-section heterogeneity and identify modest, model-dependent relationships with particular vascular and stromal surroundings. Independent-heart replication is still needed.

## Methods and scope

Only the published 21 Xenium maturation targets form the CM outcome. Counts are normalized as log1p(10,000 times target count divided by all 403 measured-target counts). The score is the mean of 21 gene Z values using the unchanged pooled reference calibration, separately for deposited and core RNA. This preserves the score used in our earlier reports; it is not a recovery of the authors' unavailable exact Seurat AddModuleScore control set. No additional CM outcome genes, inferred gene modules or snRNA-seq maturation scores enter this analysis.

The ventricular cohort is defined by the reviewed ventricular-CM flag, broad CM identity, expression-QC pass, at least 20 core transcripts and at least 20 pre-existing outcome-excluded counts. No published-label concordance, conservative ventricular mask or manually annotated region is required. Missing nuclear area is retained. Existing source identities and their confidence/assessability rules are unchanged. Core RNA is a confidence-supported transcript-assignment layer, not nuclear-only truth.

The full-cohort continuous model adjusts for flexible log-capture and log-area effects, missing nuclear area, available CM identity, and 200-µm local regions. Published non-CM labels on otherwise reviewed ventricular CMs are not used as CM subtype identities; those cells retain unresolved CM subtype. The anatomy sensitivity adds available image-region labels and retains unclassified cases. Within-section confidence intervals use 400-µm clustered uncertainty. Core full-cohort source comparisons are corrected together; deposited RNA and anatomy models are measurement/context sensitivities.

Raw high/low summaries use within-stage quartiles of the same 21-gene score. Strict matching requires the same 200-µm grid and CM subtype, a maximum 100-µm pair separation, capture and area differences no greater than 0.25 full-stage SD, identical missing-nucleus status, and unique cells. Neighbor identities and counts do not determine matches. Matched intervals use 1,999 resamples of 400-µm pair-anchor clusters. The all-CM continuous model is the main inferential view because the strict high/low match has limited measurement overlap. This is an exploratory reanalysis, not a prospectively registered confirmation.

All sources, raw comparisons, adjusted models and unsuccessful matched tests are provided. The number of cells does not substitute for independent hearts: one observed Xenium section is available per age. The score measures a transcriptional index, and these spatial relationships remain associations.

Source: Wang et al., Spatially guided in vivo single-cell functional genomics of postnatal heart, Nature Cardiovascular Research, DOI 10.1038/s44161-026-00861-z. https://pmc.ncbi.nlm.nih.gov/articles/PMC13561863/ Gene list and exact score reconstruction checks accompany this report.
