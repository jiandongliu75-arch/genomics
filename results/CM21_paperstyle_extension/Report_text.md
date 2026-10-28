# Cardiomyocyte maturation in its vascular and stromal surroundings

## What this analysis adds

The central question is whether ventricular cardiomyocytes with higher expression of the published 21-gene Xenium maturation panel occupy different vascular and fibroblast neighborhoods. All 131,605 previously eligible ventricular CM profiles are retained across P0, P7, P14 and P21. The positive maturation module contains only those 21 directly measured genes; other measured targets supply normalization and background controls. The existing snRNA-seq-supported annotations identify neighboring cell populations; they do not supply or impute the maturation genes.

The strongest result is at P7. Higher-scoring CMs are closer to capillary endothelial cells, have more capillary profiles nearby and have a greater fraction of their segmented boundary close to capillary profiles. At P14, nearness and boundary apposition show the same direction. These are small spatial associations, but agreement across several measurements makes the vascular relationship more credible than a single correlation between averaged neighborhood rings.

## The 21-gene pattern survives changes in individual genes and controls

The paper calculates a module score by subtracting expression-matched control genes from the mean normalized expression of the maturation genes. We reconstructed this approach with explicitly recorded settings suited to the 403-target panel and averaged 50 control draws to reduce dependence on random control selection. The authors' original controls, exact score object, anatomical masks and graph are unavailable, so this is a paper-style reconstruction rather than an exact replication.

Removing any one maturation gene leaves the cell ranking strongly similar to the full panel: the lowest rank correlation across stages and RNA views is 0.960. Each omission removes that gene from both the maturation module and the control pool. Changing the random controls also has little effect: the lowest correlation between an individual draw and the ensemble is 0.990. Alternative bin settings, including a setting permitting 100 controls per gene and a pool excluding all 21 maturation genes, preserve high rank agreement. Thus the result is not attributable to one naming gene or one favorable control draw (Fig. 1).

There is an important difference from our earlier report. The previous score standardized each maturation gene separately before averaging. Its rank correlation with the paper-style score is only 0.74–0.81. The two approaches therefore identify overlapping but different higher-scoring CM groups. With the paper-style score, the raw P7 comparison and the adjusted comparison now agree in direction. This explains why the earlier raw P7 result differed; it does not recover the authors' exact score or establish that one numerical index is functional maturity.

## The capillary relationship is clearest at P7 and supported at P14

The direct high-versus-low comparison uses score quartiles within each stage. At P7, the higher-scoring group is closer to capillary profiles on average. P14 shows the same direction. These are descriptive means across the section, before matching or local adjustment.

| Stage | Lower-score quartile: mean capillary distance, µm | Higher-score quartile: mean capillary distance, µm |
|---|---:|---:|
| P0 | 58.5 | 46.1 |
| P7 | 48.7 | 42.4 |
| P14 | 41.5 | 39.3 |
| P21 | 54.7 | 52.6 |

The locally adjusted analysis compares CMs within 200-µm regions and the same available CM subtype, while accounting for RNA capture, cell area and missing nuclear area. A positive coefficient means a higher score with greater nearness, more capillary profiles, or more boundary apposition. The three measurements are related descriptions of tissue geometry; they do not represent independent mechanisms.

| Stage | Capillary nearness | Capillary count within 50 µm | Boundary apposition within 2 µm |
|---|---:|---:|---:|
| P0 | +0.030; q=0.125 | +0.032; q=0.0905 | -0.013; q=0.211 |
| P7 | +0.025; q=0.0012 | +0.019; q=0.0225 | +0.011; q=0.0358 |
| P14 | +0.020; q=0.0225 | +0.011; q=0.145 | +0.014; q=0.0225 |
| P21 | +0.012; q=0.0905 | -0.002; q=0.723 | +0.009; q=0.125 |

Coefficients are standardized associations, and q values correct the 12 whole-cell vascular comparisons together. At P7, all three pass that correction. At P14, nearness and apposition pass; the count comparison does not. The P0 and P21 comparisons do not pass the same correction. Different statistical support across ages does not demonstrate when a biological relationship begins or ends (Fig. 2).

The capillary-nearness coefficient remains positive after each individual maturation gene is omitted, with every tested control setting, and in boundary-supported RNA. At P7, alternative local grid definitions also retain a positive relationship with confidence intervals above zero. At P14, the direction remains positive, although uncertainty reaches zero in one grid sensitivity. Distance along the tissue-screened cell graph and distance between segmented boundaries give the same direction at P7 and P14. The graph avoids some connections across empty spaces but is not a reconstruction of perfused vessels.

Capillary annotation is a remaining source of uncertainty. At P7, the positive association persists when using the broader public capillary labels after EC QC. Requiring strict cell-and-nucleus subtype agreement reduces the source set from 1,464 to 268 profiles and weakens the estimate to +0.008 SD, with an interval spanning zero. At P14, the broader public-label estimate is weaker, and the strict 202-profile estimate is inconclusive with a slightly negative point estimate. The evidence is therefore strongest for the reviewed and broader P7 capillary maps; it is not invariant to every source definition. Sparse, highly selected source sets are not independent ground truth for the whole network.

## Nearness is more informative than simply counting nearby cells

More capillary profiles can reflect a generally denser tissue neighborhood. After including total nearby resolved cell count, the P7 capillary-count interval includes zero, whereas capillary nearness retains a positive relationship. The P14 nearness relationship also remains positive in this sensitivity. The present evidence therefore favors a relationship with local vascular arrangement over a claim that any increase in nearby EC number predicts a higher maturation score.

Boundary apposition provides a complementary observation: a higher score accompanies a greater fraction of the CM outline lying within 2 µm of capillary profiles at P7 and P14. These are native segmented boundaries sampled at 64 evenly spaced positions. They are not confirmed membrane contacts, and the boundaries themselves can be inaccurate. The boundary-supported RNA sensitivity challenges transcript ownership while keeping the same native geometry; it does not independently validate or replace the segmentation.

## The neighboring subtype matters, but fibroblast effects remain less settled

Capillary ECs and endocardial-like ECs do not show the same relationship. In the paper-style whole-cell analysis, endocardial-like EC nearness accompanies lower scores at P7, P14 and P21, while capillary nearness accompanies higher scores at P7 and P14. This distinction argues against treating all endothelial cells as equivalent. Endocardial-like ECs also mark different anatomical environments, so regional position remains a plausible explanation. These observations do not mean that endocardial cells inhibit maturation.

Fibroblast relationships are more dependent on measurement and score definition. Lsamp-associated fibroblast nearness accompanies lower scores at P14 and P21 in whole-cell and boundary-supported RNA, consistent with the direction in the preceding 21-gene analysis. This negative direction also persists with the cell-supported, boundary-supported and crosschecked fibroblast assignments. Col8a1-associated proximity is negatively associated in the new whole-cell score screen, but the boundary-supported RNA intervals include zero. It should therefore remain a boundary- or measurement-sensitive finding, rather than a conclusion that Col8a1 fibroblasts delay maturation. Postn and Gfpt2/Pi16 nearness do not pass the new whole-cell source-screen correction where assessable. In particular, the earlier modest P7 Gfpt2/Pi16 result is not consistently recovered across scoring choices.

The full source screen also includes pericytes, smooth muscle, immune and epicardial cells, and supported EC subtypes. Results and score/segmentation disagreements are provided together (Fig. 3). Source categories overlap with their parent lineages and are not independent biological explanations. Unknown or insufficiently supported source families remain not assessable; they are not counted as absent.

## Does fibroblast context change the capillary relationship?

We tested whether the association between capillary nearness and the 21-gene score differs as the number of nearby fibroblasts of a particular family increases. These models include capillary nearness, fibroblast context, their interaction and total local resolved cell count, alongside the technical and local-region adjustments.

None of the 12 assessable whole-cell interactions passes correction for multiple comparisons. P7 Postn context gives a positive nominal interaction, but its corrected q is 0.162. We therefore cannot conclude that a Postn-rich, Gfpt2/Pi16-rich or Col8a1-rich environment strengthens or weakens the capillary relationship. This is a direct test of the proposed combined vascular–matrix niche, and its present answer is inconclusive rather than positive.

## Biological interpretation and next experiment

The most defensible interpretation is that CM transcriptional maturation and local coronary organization are associated, with the clearest evidence in the P7 section and more source-definition-sensitive support at P14. Nearness and operational boundary apposition are more consistent than capillary abundance alone. The score is heterogeneous within a stage, and this heterogeneity relates to the identity and arrangement of neighboring cells. Capillary-label uncertainty prevents calling the vascular relationship fully resolved.

These findings make local vascular integration a reasonable hypothesis to test. They do not distinguish blood flow, oxygen delivery, contact-dependent signals, matrix organization or shared regional development. A useful next experiment would compare matched ventricular regions in independent P7 and P14 hearts, measuring the same 21-gene score alongside perfused capillaries, validated CM–EC boundaries and a functional CM measurement. Earlier P3/P5 and intermediate P10 samples would establish how these relationships develop; they cannot be supplied by interpolation from the existing four stages.

## Methods and scope

Counts are normalized as log1p(10,000 times target count divided by all 403 measured-target counts). Only the 21 published maturation genes enter the positive module. The primary paper-style score uses eight rank-balanced expression bins, 20 controls per module gene sampled without replacement, the unique control union within each draw, and the mean of 50 independently seeded draw scores. Controls may include module genes, as in the default Seurat pool. Bins use frozen pooled all-profile expression means and are fixed across stages within each RNA view. Rank binning is explicit and is not identical to Seurat's implementation. Alternative settings use three bins/100 controls, 24 bins/five controls, and eight bins/20 controls excluding the 21 panel genes. These settings challenge the reconstruction; none is asserted to recover the authors' original pool or bins. The 403-target panel alone cannot support 100 distinct controls in each of 24 balanced bins, which is a feasibility issue for that parameter combination, not evidence that the authors made an error.

The primary view is deposited whole-cell Xenium RNA, matching the paper's measurement context. Confidence-supported core RNA is a transcript-assignment sensitivity and is not nuclear-only truth. The prior all-QC ventricular cohort is unchanged: no anatomical exclusion, public-label agreement or visible-nucleus requirement is added. Non-CM source labels and their eligibility/assessability rules are inherited from the reviewed reference. Nuclear-supported EC labels and nuclear-optional multi-gene fibroblast proposals remain separate evidence layers. No new segmentation or new subtype clustering was performed here.

Local regressions use view-specific non-panel molecule counts and cell area as cubic splines, missing nuclear area, available CM subtype and 200-µm grid-by-CM-state fixed effects. Unknown CM subtype and unknown candidate anatomy are retained. Confidence intervals use 400-µm clustered uncertainty with an absorbed-group degrees-of-freedom correction. Primary vascular q values use 9,999 null-imposed clustered-score sign draws, followed by BH across 12 tests; this is not vessel relocation. Non-CM nearness and fibroblast-context tests have separate BH families of 48 and 12. All uncertainty concerns these observed sections, not independent animals. Supporting tests include alternative score settings, all 21 gene omissions, core RNA, existing candidate anatomy, local grid changes, 800-µm uncertainty blocks, 30-µm counts and alternative graph/boundary distances. Supporting comparisons are not additional discoveries.

Source-label sensitivities hold the CM score and cohort fixed while substituting EC-QC-supported public capillary labels, strict cell-and-nucleus EC calls, or cell-supported/boundary-supported/crosschecked fibroblast family calls. Their estimates and source coverage are provided without selecting a favorable source definition. The native boundaries and tissue-screened graphs are inherited, not newly expert-validated anatomical masks. Independent validation re-reads count matrices, reconstructs all 50 control draws, verifies complete cohort retention, recomputes local capillary regressions and rereads original boundary vertices for sampled CMs. One observed section per stage prevents animal-level replication, and no result establishes a causal effect or functional CM maturation. Source: Wang et al., Nature Cardiovascular Research, DOI 10.1038/s44161-026-00861-z; https://pmc.ncbi.nlm.nih.gov/articles/PMC13561863/.
