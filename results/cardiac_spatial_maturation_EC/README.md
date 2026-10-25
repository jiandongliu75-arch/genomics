# Cardiac spatial organization and ventricular maturation near active ECs

**The P0 and P7 sections do not provide convincing evidence that ventricular cardiomyocytes near proliferating endothelial cells have lower overall 21-gene maturation scores.** The spatial distributions of CM, fibroblast and EC subtypes are more clearly differentiated than this proposed local maturation relationship.

[Ten-page report](cardiac_spatial_maturation_EC_report.pdf) · [Association figure](EC_activity_CM_maturation_associations.png) · [Near-cycling-EC contrasts](near_cycling_EC_maturation_contrasts.tsv) · [All models](association_coefficients.tsv) · [Supported-subtype neighborhoods](subtype_lineage_neighborhoods_30um.tsv) · [Validation](validation.json)

**Methods review:** [what the models estimate, how the FB/EC pipelines differ, and a small model-dependent P0 signal](METHODS_REVIEW.md). The primary results are unchanged; this review makes adjustment and cohort limitations explicit.

## Spatial organization

The maps preserve all native positions and the preceding broad identities. They display the three requested lineages—cardiomyocytes, fibroblasts and endothelial cells—over the remaining cardiac profiles. Gray fine identities are retained. Strict ventricular CM selection covers 149,623 profiles; broad CM labels cover 183,222. Fibroblast family maps use the preceding strict A+B cell/nuclear agreement. EC maps explicitly use nuclear-supported labels in the eligible 29,363-EC cohort; those labels are distinct from the stricter cell+nuclear subtype view, which is also summarized in the tables.

| Population | Spatial finding supported by the current maps |
|---|---|
| Working/cycling ventricular CM candidates | Distributed in myocardial territories, rather than occupying an exclusively proliferating-EC-associated compartment. Only supported subsets of the fine RNA groups can be placed. P0 cycling and neonatal-working candidates have similar broad30µm neighbor compositions. |
| Slit2-associated conduction-program CM candidates | Closer to endocardial proxies: median distances about 14.5,17.5 and17.6µm at P0/P7/P14 for the prioritized candidates within strict ventricular CMs. This supports an inner/lining-associated position, not verified Purkinje identity. |
| Slit2/Mpped2 regional-electrical CM candidates | Selected candidates at P14/P21 have higher21-gene panel scores, with median 3.31/3.32. Fine labels use overlapping expression information, so these differences are descriptive rather than independent subtype validation. |
| Ankrd1-associated CM candidates | Localized remodeling candidates at P14, with median candidate-arterial-envelope distance188µm in this prioritized subset. They are not confined to an exclusively peri-arterial compartment. This descriptive subset does not replace the earlier whole-cohort vessel models. |
| Postn-associated fibroblasts | Early, widespread CM-containing/interstitial and stromal territories. P0/P7 strict-family profiles have approximately 58%/46% CM neighbors at 30µm. Developmental remodeling is not pathological fibrosis. |
| Gfpt2-associated fibroblasts | CM/EC-containing interstitial territories; relatively fewer fibroblast neighbors than Postn at P7 and fewer smooth-muscle neighbors than Col8a1 at P14. |
| Col8a1-associated fibroblasts | CM/EC-containing interstitium with a partial vascular-associated component, clearest at P14. The preceding geometry-based comparison places18.8% within 50µm of candidate arterial envelopes versus7.0% of Gfpt2 profiles; most remain outside that envelope. |
| Lsamp-associated / Kirrel3-associated fibroblasts | Too few strict spatial calls to establish reliable separate niches. Do not assign outflow tract or valve anatomy from a marker or permissive proposal. Separate Col8a1/Svep1 versus Col8a1/Hmcn2 and Lsamp/Slit3 versusLsamp/Col4a5 RNA substates remain unresolved spatially. |
| Capillary-like ECs | Interspersed with CMs through myocardial territories. Strict-type profiles have approximately55–76% CM neighbors at 30µm across the sections. |
| Endocardial-like ECs | More concentrated in lining-associated territories, with relatively more EC/stromal neighbors than capillary-like profiles. Manually traced chamber anatomy is unavailable. |
| Arterial-like ECs | Associated with smooth-muscle-containing vascular structures. Smooth-muscle neighbors comprise approximately24%,38% and28% at P0/P14/P21, versus approximately0.1–0.9% around capillary-like ECs. P7 conservative arterial typing remains unresolved. |
| Venous-like / lymphatic-like EC candidates | Limited local candidate assignments, especially at P14. Sparse lymphatic calls and incomplete panel coverage do not support an exhaustive venous/lymphatic network or unbiased absence claims. |

The named fibroblast observations incorporate the [preceding confidence-qualified physical-neighbor analysis](../fibroblast_names_and_neighborhoods/README.md), using the same frozen broad atlas and family decisions. The new table directly recalculates30µm lineage neighborhoods for CM candidates, strict fibroblasts, and strict/nuclear-supported EC types. Sparse groups are marked rather than interpreted as distinct anatomical territories. The EC strict and nuclear-supported views overlap and must not be added together.

## Ventricular maturation and the P0/P7 hypothesis

The outcome is the **mean log1p(CP10k) of the exact 21 measured Figure4b genes**, not a newly trained maturity classifier or a claim of reproducing the authors' precomputed index. The inherited score increases across the four sections:

| Stage | Strict ventricular CM profiles | Mean 21-gene score |
|---|---:|---:|
| P0 | 19,304 | 2.532 |
| P7 | 39,383 | 2.675 |
| P14 | 42,796 | 2.923 |
| P21 | 48,140 | 2.983 |

Fixed expression-matched-control scores also increase. These age comparisons remain descriptive because there is one spatial section per age and capture/section geometry differ.

An EC is **nuclear cycling-marker-positive** when at least two distinct genes among Mki67,Top2a,Cenpf have nuclear molecules. This is marker support, not completed division. Among the selected eligible EC cohort, counts are763/3,732 at P0,1,185/8,944 at P7,114/8,086 at P14 and91/8,601 at P21. These confidence/capture-selected Xenium fractions are not directly comparable with snRNA fractions or whole-organ proliferation rates.

Primary analyses require at least3 eligible ECs within 50µm and at least5 resolved cardiac neighbors, retaining CMs with missing nuclear-area measurements. The near-cycling comparison uses CM centroids within 30µm of a nuclear cycling-positive EC; comparison CMs are other eligible CMs in the same section, rather than another age or region selected for a desired maturity result.

| Stage | Near CMs | Comparison CMs | Adjusted near−comparison maturation | 95% within-section spatial-block interval |
|---|---:|---:|---:|---:|
| P0 | 7,012 | 10,731 | +0.0056 | −0.0073 to +0.0187 |
| P7 | 11,042 | 25,203 | +0.0062 | −0.0009 to +0.0125 |
| P14 | 730 | 27,214 | +0.0035 | −0.0213 to +0.0293 |
| P21 | 445 | 25,677 | −0.0375 | −0.0659 to −0.0029 |

Differences are in the original mean log-expression score units. P0/P7 near/far differences are tiny, with intervals including zero. The primary local continuous nuclear cycling module, local cycling-positive EC fraction, inferred proliferation module and relative proximity to cycling versus marker-negative ECs give the same broad early-stage conclusion. Removing CM cycling/Slit2 from adjustment reveals a small negative P0 coefficient for continuous nuclear EC cycling (−0.01686SD; block interval−0.03294 to−0.00104), while P7 and the binary proximity comparison remain uncertain. This model dependence prevents claiming that every sensitivity is null; it does not establish a robust P0/P7 immature-CM niche. See the [methods review](METHODS_REVIEW.md).

**Across-age co-trending is different from within-age spatial association.** Younger sections show lower CM panel scores and more EC cycling-marker support, but that does not establish that individual CMs closer to active ECs are less mature withinP0/P7.

## Sprouting inference and important qualifications

The six-gene inferred sprouting signature has a tiny positive, rather than negative, association with the native CM maturation score at P0/P7: approximately+0.016/+0.019 maturationSD per local-exposureSD. These are small descriptive effects, not verified biological sprouting niches. Their added buffered spatial predictiveR² is at most 0.00043 at these ages; early proliferation exposures also add very little predictive information beyond the baseline.

OnlyAngpt2 from the original sprouting signature is measured by Xenium. The missing-five component is an explicit diagnostic. A new sensitivity reprojects the RNA programs after withholding **all 21 CM maturation genes plus all sprouting/cycling targets from EC alignment features and their normalization**. This leaves 353 shared predictors. Early RNA opposite-library recovery of sprouting remains weak, especially P7. The sensitivity avoids direct reuse of the maturity outcome as an EC predictor; it does not independently validate missing expression or erase inherited outcome-dependent broad cell selection. Actual Xenium count matrices are untouched.

**Exploratory component finding:** the P0 sarcomere-only panel component has a modest negative coefficient with local nuclear cycling (approximately−0.045SD; within-section interval−0.060 to−0.024). The overall 21-gene panel and P7 do not show that pattern, and the near-versus-other binary comparison does not show a P0 sarcomere difference. These four component tests are exploratory and lack animal replication; they are not evidence of global CM immaturity near cycling ECs.

**P21 qualification:** the small negative native-count proximity contrast attenuates using CM boundary-core counts and becomes uncertain in the author-compatible ventricular subset. It is boundary-sensitive and should not be promoted to a reliable late-stage niche. A positive or negative inferred score association is not causal evidence that ECs regulate maturation.

## Covariates, sensitivities and uncertainty

Each exposure is fitted separately, within each section. The primary conditional model includes CM RNA depth excluding the outcome panel and previous Ankrd1-related outcomes, native area, missing nuclear area, local resolved-profile density, eligible EC number, nearest eligible EC distance, FB/EC/mural/immune proportions, local EC nuclear capture and CM-marker fraction, EC subtype proportions, CM cycling and Slit2 scores, XY splines/interaction and existing arterial-envelope/endocardial-proxy distances. These are incomplete tissue-structure proxies rather than manually established chamber/transmural anatomy.

Sensitivities cover30/100µm neighborhoods;200/800µm blocks; a model without spatial/anatomical terms; removal of CM state covariates; fixed expression-matched outcome controls; removal ofMyl2 from the outcome; resegmented/boundary-core CM outcomes; native/core EC cycling scores; gene-balanced sprouting; author-ventricular compatibility; and a new EC projection excluding all maturity/sprouting targets. Rebuilding neighborhoods using only individual EC profiles with nuclear CM-marker fraction≤5% tests a more restricted EC source; it is an operational restriction, not proven RNA purity. The proposed≤2% neighborhood-average restriction yields fewer than 200CMs per age and is explicitly recorded as unfit rather than silently reported as a negative result.

Bootstrap intervals use 256 resamples of 400µm spatial blocks for primary models, with128/96 draws for sensitivity models. They describe within-section sampling variability, not independent animals. Buffered five-fold prediction uses 400µm blocks and removes training CMs within 100µm of held-out CMs; overlapping EC neighborhoods and a finite set of shared EC sources still limit independence. Coefficients use standardized exposure and outcome units; binary near-versus-other contrasts are additionally exported in original outcome units. Many exploratory models are reported transparently; no multiplicity-adjusted animal-level significance is asserted.

## Maps and downloadable data

| Stage | Six-panel spatial overview | CM maturation | EC activity | CM/EC 50µm neighborhoods |
|---|---|---|---|---|
| P0 | [PNG](P0_cardiac_spatial_overview.png) | [TSV.gz](P0/ventricular_CM_maturation.tsv.gz) | [TSV.gz](P0/EC_activity_annotations.tsv.gz) | [TSV.gz](P0/CM_EC_neighborhoods_50um.tsv.gz) |
| P7 | [PNG](P7_cardiac_spatial_overview.png) | [TSV.gz](P7/ventricular_CM_maturation.tsv.gz) | [TSV.gz](P7/EC_activity_annotations.tsv.gz) | [TSV.gz](P7/CM_EC_neighborhoods_50um.tsv.gz) |
| P14 | [PNG](P14_cardiac_spatial_overview.png) | [TSV.gz](P14/ventricular_CM_maturation.tsv.gz) | [TSV.gz](P14/EC_activity_annotations.tsv.gz) | [TSV.gz](P14/CM_EC_neighborhoods_50um.tsv.gz) |
| P21 | [PNG](P21_cardiac_spatial_overview.png) | [TSV.gz](P21/ventricular_CM_maturation.tsv.gz) | [TSV.gz](P21/EC_activity_annotations.tsv.gz) | [TSV.gz](P21/CM_EC_neighborhoods_50um.tsv.gz) |

Use `stage_cell_id` to combine ages; native cell IDs can repeat between sections. Supporting tables include [fine CM candidate summaries](fine_CM_candidate_spatial_summary.tsv), [EC restriction coverage](EC_CM_marker_restriction_coverage.tsv), [outcome-withheld RNA benchmarks](CM_panel_withheld_RNA_library_benchmarks.tsv), [spatial prediction](spatial_blocked_prediction.tsv), [coverage](neighborhood_coverage.tsv), and [maturation stage summaries](maturation_stage_summary.tsv). Per-stageNPZ files preserve the new withheld-alignment genes, actual reference target expressions and every RNA-neighbor index/weight.

Reproduce with `OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 /workspace/heart-atlas/venv/bin/python analysis/cardiac_spatial_maturation_EC.py --part all`, then `analysis/report_cardiac_spatial_maturation_EC.py` and `analysis/validate_cardiac_spatial_maturation_EC.py` in the same environment. Validation independently checks marker arithmetic, source hashes, unchanged CM scores and coordinates,256 physical-neighbor examples across all three radii, new module projections reconstructed from weights, primary coefficients reconstructed by partial regression and figure/PDF integrity. These are implementation checks, not independent anatomy or biological replication.

The deposited Xenium sections and source images/counts are attributed to the paper's authors, as in the [source atlas](../whole_heart_spatial_atlas/README.md). No original expression matrix, cell identity or physical coordinate is replaced by this analysis.
