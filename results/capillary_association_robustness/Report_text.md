# Strengthening the capillary–cardiomyocyte association

## Main conclusion

The additional tests strengthen a modest association between capillary proximity and CM transcriptional maturation within the P7 section. The result is similar when distance is measured from cell centers, from the sampled CM boundary, or along a tissue-screened graph. It is also stable to changes in the size and position of the local comparison regions. These checks reduce concerns that the finding depends on one geometric definition or a particular grid boundary.

They do not establish that the association is strong, that every CM subtype behaves identically, or that capillaries cause maturation. A tightly matched near–far CM comparison remains inconclusive. Independent hearts are still the most important missing evidence.

## Does the result depend on how distance is measured?

The primary comparison uses 30,485 ventricular CMs, or 30,484 with a reachable graph path. Cells are compared within the same CM subtype and 200-µm region, with flexible adjustment for RNA capture and cell size, missing nuclear area, morphology and cavity/surface distances. All three definitions give similar small positive associations:

| Distance definition | Score SD per SD of nearness | 95% spatial-cluster interval |
|---|---:|---|
| Cell-center distance | +0.021 | +0.012 to +0.030 |
| Sampled CM boundary distance | +0.019 | +0.012 to +0.026 |
| Tissue-screened graph path | +0.020 | +0.012 to +0.029 |

The boundary analysis helps address a specific concern: larger CMs can have distant centers even when their surface is close to an EC. Its agreement with the centroid analysis makes this geometric explanation less likely to account for the whole association. Boundary distance is still approximate: it uses 64 sampled points on the deposited CM polygon and capillary territories, rather than independently validated membrane positions. The graph-path analysis avoids some straight-line shortcuts through cavities, but is not a true tissue geodesic or a measure of vascular flow.

All three associations also pass the prespecified cluster-score sign-flip sensitivity: none of 9,999 simulated sign-flip statistics was as extreme, giving a Monte Carlo P of 0.0001 and the same three-test BH q. This is a conditional statistical check using 400-µm spatial clusters, not random relocation of vessels or replication across hearts. Residual correlations across clusters and unmeasured anatomy remain possible.

## Could one region, one subtype or local tissue density explain it?

For each distance definition, I repeated the comparison using 100-, 200- and 400-µm regions and four grid positions at each size. All 36 estimates and their clustered intervals remain positive. This reduces concern that the initial result was created by the chosen 200-µm grid. Removing any single 400-µm region from the fitted partial-regression sums also leaves a positive estimate. This last check is an influence diagnostic with nuisance residuals held fixed, not independent held-out validation.

The association persists when capillary sources are defined from broad EC support plus public capillary labels without requiring nuclear subtype agreement: +0.024 SD, with an interval of +0.016 to +0.031. It therefore does not rely entirely on the stringent nuclear-supported capillary subset.

Additional conditional models include local cell density, capillary-source capture and capillary counts, or distances to fibroblasts, noncapillary ECs and pericytes. The capillary-nearness coefficient remains positive. With local density, source capture and capillary count held similar, it is +0.017 SD (+0.005 to +0.028). These comparisons make a purely generic relationship with cellular crowding less convincing. They do not show that capillaries act independently of the rest of the vascular or stromal environment; some of these controls can represent parts of the same biological relationship.

Subtype-specific point estimates are positive for working ventricular, Ankrd1-associated, Ddc-associated, cycling and Slit2-associated CMs. Several subtype intervals cross zero, especially in smaller or spatially restricted groups. A significant result in one subtype and a nonsignificant result in another does not establish a difference between subtypes. These analyses are secondary and do not justify a claim that the relationship is universal or stronger in one particular subtype.

## What does the stricter local match show?

I compared CMs within 25 µm of a QC-supported capillary center with CMs at least 50 µm away. A pair had to share CM subtype, morphology and 200-µm region, lie within 100 µm of each other, and have closely similar capture and cell area. Matching did not use the maturation score, and no CM was reused. This yielded 916 pairs across 62 spatial clusters.

The mean 21-gene difference is +0.021 SD, with a spatial-block interval of -0.040 to +0.079. Deposited RNA and equal-molecule scores also have positive estimates with intervals crossing zero. This comparison is compatible with the small regression estimate, but it cannot independently confirm a near–far difference. It uses a smaller selected subset and a specific distance contrast, so its wider uncertainty should be retained rather than explained away or hidden.

## What would strengthen the biological conclusion most?

First, repeat the analysis in independent P7 hearts and additional early stages using the same frozen score and distance definitions. The heart, rather than the number of cells or image fields, should be the biological replication unit. This is the most important improvement and cannot be supplied by further resampling of this section.

Second, annotate compact and trabecular myocardium, ventricular chamber, papillary regions and large-vessel surroundings independently from morphology, preferably without viewing maturation scores. Our current local comparisons and morphology proxies reduce regional confounding, but do not recover the original authors' exact anatomical masks. Prespecified analyses within these anatomical regions would make a regional developmental explanation easier to distinguish from a more local vascular relationship.

Third, validate the association against protein and CM function. Sarcomere organization, calcium handling and contraction would test whether the transcriptional association accompanies functional development. Perfusion tracing and direct oxygen measurements would characterize vascular support. PTPRM localization and cell-specific perturbations could then test a candidate mechanism. Neither mechanism should be inferred from the distance association itself.

The strengthened statement is: At P7, capillary proximity has a small association with ventricular CM transcriptional maturation that is robust to several distance definitions, local comparison scales and source-confidence rules. A tightly matched near–far comparison remains inconclusive, and independent-heart replication is needed.

## Methods and reproducibility

This follow-up retains the prior P7 strict ventricular CM cohort, 21-gene frozen-reference Xenium score, corrected RNA views and capillary definitions. No new segmentation, imputation, cell taxonomy, samples or causal inference was introduced. Core RNA is a confidence-supported assignment layer, not nuclear-only truth; anucleate CM sections remain eligible. The three primary geometry tests were chosen before examining their results. Other source definitions, grid, subtype, density and matching analyses are reported as sensitivities, without selecting favorable results as additional primary tests.

The score and exact original author masks/graph remain unavailable as author objects; our score is a directly measured gene-Z mean rather than an exact recreation of AddModuleScore. Regression uses flexible capture/area splines, within-region and CM-type centering, and 400-µm clustered uncertainty. Sign-flip inference uses null-imposed clustered partial-regression scores and assumes approximate cluster independence and symmetry; it does not preserve or randomize vessel topology. Cell-center and boundary distances are two-dimensional measurements. Graph distance sums centroid edge lengths on the previously tissue-screened Delaunay graph. Capillary count and general tissue controls define narrower conditional questions and are kept separate from the primary proximity estimate.

The local match uses fixed thresholds of at most 25 µm versus at least 50 µm from capillary centers, a 100-µm maximum pair separation, the same 200-µm grid/CM state/morphology, and at most 0.25 full-cohort SD differences in log capture and log area. Pair intervals use 1,999 resamples of 400-µm spatial clusters. They describe within-section uncertainty and do not substitute for independently sampled animals.
