# Capillary proximity and cardiomyocyte development: a targeted retest

## What the retest shows

CMs closer to capillaries tend to have slightly higher 21-gene maturation scores when we compare cells of similar identity, anatomy and RNA capture. At P7, the association is +0.051 (+0.040 to +0.061) SD of score per SD of nearness, using the QC-supported capillary definition. The relationship becomes smaller when we compare CMs of the same subtype within 200-µm regions and allow nonlinear capture and cell-size effects: +0.021 (+0.012 to +0.030). It persists in deposited RNA, after equal-molecule subsampling and after leaving individual maturation genes out. These checks support a modest local transcriptional relationship with vascular proximity.

The raw P7 relationship is not consistently positive. The association appears after capture and tissue-context adjustment, so it should not be presented as a strong, unconditional rule that every CM near a capillary is more mature. Many nearby CMs still have different developmental states. We also cannot recover the authors' exact score-bearing object, anatomical masks or original graph. This is a new test of the biological question, rather than an exact numerical reproduction of Figure 4e.

| Age | Ventricular CMs | Adjusted association, SD (95% interval) |
|---|---:|---|
| P0 | 17,272 | +0.054 (+0.035 to +0.073) |
| P7 | 30,485 | +0.051 (+0.040 to +0.061) |
| P14 | 33,659 | +0.037 (+0.029 to +0.045) |
| P21 | 37,949 | +0.014 (+0.005 to +0.022) |

All four QC-source associations pass the eight-test primary within-section correction, which also includes the public capillary labels. The intervals describe spatial variation in one observed section per age. They do not establish reproducibility across animals or a developmental increase or decrease in effect.

## Why the paper's relationship can look strong

The paper summarizes CM scores by graph-neighborhood level. In our P7 reconstruction, the correlation between mean score and neighborhood level is approximately -0.90, while the correlation across individual CMs is only about -0.044. The means show an organized spatial trend, but the individual-cell relationship is weak. Averaging removes much of the cell-to-cell variation; ten neighborhood means are not ten independent hearts. Farther rings are also sparse: the reconstructed tenth ring contains only two strict ventricular CMs, compared with 14,160 in the first. The trend should therefore be assessed with cell numbers and uncertainty, rather than judged from the correlation of means alone. Euclidean distance, graph steps and boundary apposition describe different aspects of vascular integration and need not give the same coefficient.

The new continuous apposition regression is positive at P7, but our earlier tightly matched boundary comparison remains compatible with little or no difference in overall maturation score. The same matched-cell contrasts were independently recomputed here. Thus, vascular proximity is better supported than a precise claim that physical contact itself drives maturation. Matching additionally held local capillary area, cell shape and capture tightly similar and selected a smaller CM subset; that comparison asks a narrower question.

## Does PTPRM offer a mechanism?

PTPRM is a reasonable contact-associated candidate. The atlas shows PTPRM protein at EC–CM boundaries and reports increased cell-cycle RNA and reduced Meis1 and maturation-related RNA after Ptprm perturbation. This is more informative than a spatial correlation alone. Because PTPRM is expressed by both CMs and capillary ECs and participates in homophilic interactions, however, the experiment does not by itself distinguish EC-to-CM signaling from CM-intrinsic or CM-to-CM functions.

Our corrected snRNA-seq reference supports Ptprm expression in CMs and capillary ECs in both libraries at every age. In P7 Xenium, CM Ptprm expression is associated with higher overall and metabolic-regulatory scores and lower cycling RNA. The relationship is not uniform across maturation components: sarcomere RNA has the opposite direction. With flexible capture controls, the overall-score coefficient for CM Ptprm is +0.038 (+0.030 to +0.045).

The more discriminating tests do not establish a capillary source. For P7 CMs sharing a refined territory boundary with a QC-supported capillary, the coefficient for capillary Ptprm is +0.009 (-0.008 to +0.025), and that for the capillary-by-CM Ptprm expression product is -0.005 (-0.013 to +0.003). Both intervals cross zero. Neither term gives a consistent positive relationship across deposited and core RNA. This leaves PTPRM as a biologically plausible candidate, without identifying endothelial PTPRM as the explanation for the proximity association.

This distinction has a useful biological implication. In the P7 RNA reference, Ptprm is detected in at least 98.8% of capillary nuclei in each library. If endothelial PTPRM is broadly available, differences in its RNA abundance may not determine signaling strength. CM competence, protein localization at the interface or the stability of the interface could matter more. A weak capillary-RNA coefficient therefore does not rule out a requirement for endothelial PTPRM. Separate EC and CM perturbations are more discriminating than demanding a stronger spatial expression correlation.

Detection and geometry limit this comparison. Only 3,684 of 30,485 P7 CMs have a QC-supported capillary neighbor on the refined contact graph; 1,349 have detected core Ptprm on both sides. Absence of detected RNA is not absence of protein or signaling. Shared segmented territory boundaries are operational contact candidates, rather than verified membrane interactions. Meis1 is absent from the Xenium panel, so the paper's proposed downstream link cannot be directly followed spatially here.

## Could oxygen delivery explain the relationship?

Yes, but capillary proximity is an incomplete measure of oxygen supply. Supply depends on whether vessels are perfused, local flow and oxygen content, diffusion distance and the CM's oxygen demand. A capillary in a two-dimensional section can be close to a CM without adequately supplying it. Conversely, a nearby vessel can lie outside the plane of the section. More developed CMs may also recruit or retain vessels, so reverse causation remains possible.

Oxygen-dependent metabolic adaptation and PTPRM-associated contact are compatible mechanisms. Vascular growth could first establish access to the tissue, while subsequent perfusion and stable cell interfaces support different aspects of CM development. Oxygen-driven mitochondrial and redox changes could contribute to cell-cycle withdrawal, while contact-associated signals affect how CMs respond to that environment. Our data cannot assign the proximity relationship to either route.

The Xenium panel measures only Vegfa among the 16 genes in our hypoxia-response panel, and only Cox5a among the 18 genes in the broader oxidative-metabolism panel. Vegfa is regulated by several developmental inputs; it is not an oxygen meter. Nuclear RNA provides broader response programs, but those nuclei do not carry independently measured spatial oxygen exposure. Their hypoxia scores are also sparse and are standardized within library, so they cannot estimate a whole-heart change in oxygen across ages. Projecting those programs onto Xenium would produce a hypothesis map, not a measured oxygen map. We therefore have not labeled regions as hypoxic or inferred perfusion from maturation genes.

## The most useful next experiments

First establish whether the association follows vascular access, actual cell interfaces, or both, in independent P3/P5 and P7 hearts. Combine CM boundaries and PTPRM protein with endothelial staining and a perfusion tracer. A tracer can identify vessels reached by the circulation, although it does not quantify local oxygen or flow rate. Where feasible, add a direct tissue-oxygen method; phosphorescence-lifetime measurements can measure local oxygen, whereas pimonidazole mainly identifies sufficiently low oxygen and cannot resolve the full physiological gradient. Analyze comparable ventricular regions and CM subtypes, and use the heart as the replication unit.

Then separate contact from oxygen under controlled conditions. Compare CM–EC contact with a noncontact condition while measuring and matching local oxygen exposure. Perturb Ptprm separately in ECs and CMs, with appropriate controls and rescue. An EC-specific requirement under matched oxygen would favor an EC contact-associated mechanism. A stronger dependence on oxygen at fixed contact would favor oxygen support. A combined effect would support cooperation between these routes. Any perturbation that changes vascular structure or perfusion needs those changes measured before it can be interpreted as a direct contact effect.

Use several CM outcomes: cycling, metabolic adaptation, sarcomere organization, calcium handling and contraction. Cell-cycle RNA alone cannot establish completed cell division or completed maturation. At present, the defensible conclusion is that local vascular integration accompanies CM transcriptional development, with a small residual proximity association; the responsible mechanism remains unresolved.

## Methods and reproducibility

The primary retest uses frozen ventricular identities and the directly measured Xenium 21-gene panel at P0, P7, P14 and P21. The score is the existing pooled-reference mean of gene Z values after log1p CP10k normalization, not the authors' unavailable AddModuleScore object. Deposited and core transcript assignments are analyzed separately; core RNA is not nuclear-only truth. Previously valid anucleate CM sections remain eligible. No new segmentation, imputation or lineage inference was performed.

Nearness is minus log1p of the nearest capillary centroid distance; a positive coefficient means a higher score nearer capillaries. Primary adjustment includes outcome-excluded CM capture, cell area, missing nuclear area, CM identity, broad morphology and distances to the outer surface and large cavities. It does not condition on CM receptors or inferred communication, so the total proximity question is tested before adding potential mediators. Sensitivities use within-200-µm-region and CM-state comparisons, cubic splines for capture and cell size, two capillary definitions, two count views, 400/800-µm spatial clusters, operational boundary apposition, all 21 single-gene omissions and fixed-50-molecule subsampling. The three-bin expression-matched score is an additional feasible sensitivity, not the recovered authors' control selection.

Contact-associated molecular models include capillary counts and capture, nearness, neighboring-CM counts and Ptprm, receiver Ptprm, capillary Ptprm and their product. Molecular normalization excludes the 21 score genes and three cycling genes from its denominator. Both main expression terms accompany the product; their coefficients are conditional associations, not measured protein binding or mediation effects. Sixty primary PTPRM tests cover four ages, five outcomes and three molecular terms in core RNA. Secondary flexible-control analyses are reported as sensitivity results. Clustered t intervals use CR1 sandwich uncertainty conditional on this observed section, with 400-µm clusters primary and 800-µm clusters as a scale check. They do not include animal, annotation or segmentation uncertainty.

Source: Wang et al., Spatially guided in vivo single-cell functional genomics of postnatal heart, Nature Cardiovascular Research, DOI 10.1038/s44161-026-00861-z; Figures 4–5, Extended Data Figure 7 and the maturation/neighborhood Methods. Public PIP_seq code was reviewed at 102d80e6ea5de44f594461c727a4b3c6e1bc3508. https://pmc.ncbi.nlm.nih.gov/articles/PMC13561863/
