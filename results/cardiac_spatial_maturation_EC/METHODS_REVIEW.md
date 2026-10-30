# Review of methods relating ventricular maturation to non-CMs

The analyses test **within-section spatial associations between a measured CM expression score and nearby non-CM identities/programs**. They do not measure a causal effect of non-CMs on CM maturation. The fibroblast and EC analyses are separate models with different eligible CM cohorts; the recent combined maps do not constitute one joint model of all non-CM states.

## What was actually measured

For each strict ventricular CM profile, the primary outcome averages log1p(CP10k) across the exact 21 measured Figure4b genes. The denominator is RNA from the 403-gene Xenium panel, not the entire transcriptome. This is a panel-expression index rather than a calibrated physiological maturity measurement. All genes receive equal weight, which gives the nine transcription-factor genes more aggregate weight than the three sarcomere genes. The four supplied categories contain 9, 3, 5 and 4 genes. Separate-category and fixed expression-matched-control scores help test whether a composite masks a component-specific pattern.

The strict ventricular selection requires a confident broad CM assignment and measured Myl2≥2 with Myl2≥Myl4. Because Myl2 is also a score gene, the source cohort is not selected independently of the outcome. Removing Myl2 from the score is useful but does not recover ventricular CMs excluded by that selection. Missing nuclear area does not automatically remove a CM profile.

## How non-CMs were connected to each CM

| Analysis | Neighborhood exposure | Important scope distinction |
|---|---|---|
| Fibroblast families | Mean same-age RNA-reference family weights among fibroblasts within 50µm;30/100µm sensitivities | Primary association models use native-count soft assignments from the reviewed refit, including uncertain proposals. Weights are similarities, not calibrated identity probabilities. Strict nuclear/cell family maps are a separate evidence view. |
| Fibroblast matrix programs | Mean measured, normalized FB-only matrix-program expression within the neighborhood | This measures present fibroblast RNA, not deposited matrix abundance, stiffness or earlier matrix exposure. Annotation markers are excluded from the tested matrix panels; an ECM-withheld mapping sensitivity is available. |
| EC activity | Mean nuclear cycling score, fraction with≥2 nuclear cycling genes, or RNA-inferred program score among eligible ECs within 50µm | Direct proliferation-marker evidence and inferred sprouting have different reliability. Only Angpt2 of the original six-gene sprouting signature is assayed in Xenium. |
| EC proximity | Within30µm of an eligible EC with≥2 nuclear cycling genes; relative distance to cycling versus marker-negative eligible ECs | Distances are between original profile centroids. These are not cell-boundary distances or proof of contact. Marker-negative ECs may include false negatives. |

The primary FB analysis requires≥3 FBs and≥5 resolved cardiac neighbors. The EC analysis requires≥3 eligible ECs and≥5 resolved neighbors. Both also require an available arterial-envelope distance. Thus the models describe selected neighborhoods, not all ventricular CMs. At50µm the FB cohorts contain 8,015/24,802/24,905/21,465 CMs and EC cohorts 17,743/36,245/27,944/26,122 at P0/P7/P14/P21. The eligible EC source comprises 29,363 of 69,617 broad native EC profiles. QC selection and sparse marker detection can change the observed exposure as well as which CMs are analyzed.

Euclidean neighborhoods may connect profiles across a lumen or tissue gap and may poorly represent adjacency for large elongated CMs. Neither possibility has been eliminated by the centroid model. Nuclear RNA improves transcript ownership but does not itself establish contact geometry.

## Statistical question and uncertainty

Each age is analyzed separately. Linear regressions estimate an exposure coefficient conditional on capture/area, local composition and density, CM cycling/Slit2 expression, and spatial/anatomical proxies. EC models additionally adjust for eligible EC capture, CM-marker signal and subtype composition. FB-family coefficients are relative composition contrasts, with a reference family omitted when weights sum to one. They are not independent changes in absolute fibroblast abundance.

Useful safeguards already implemented include withholding the 21 CM score genes from primary fine FB transfer, repeating EC projection without those genes, raw/separate-library reference checks, fixed-control and boundary-core outcomes, spatial-block inference, and buffered prediction. These reduce specific sources of circularity or technical bias without independently validating the biological model.

This full model asks whether an exposure adds an association **beyond the adjusted tissue and CM state**. If a non-CM influence acts through CM cycling or a regional developmental program, those covariates could remove part of the relationship of interest. Conversely, omitting anatomy can leave regional confounding. A conditional result and an overall local developmental association answer different questions; neither model alone establishes a total causal effect.

We also compare held-out prediction with/without the exposure using five 400µm spatial-block folds and a 100µm training buffer. A small incremental R² means little extra prediction beyond that baseline; it does not demonstrate biological irrelevance. At 100µm neighborhood radius, this 100µm buffer would not exclude every shared source-cell neighborhood, although primary prediction uses 50µm neighborhoods.

Coefficient intervals resample 400µm blocks: 96 draws in the FB analysis, 256 in primary EC models, and generally 96/128 in sensitivity models. These intervals reflect within-section variation. They do not propagate every reference/identity/segmentation uncertainty or substitute for independent animals. Only one Xenium section per age is available. Many model/category tests are exploratory, and weak intervals should not be interpreted using a sharp zero-crossing rule.

## A model-dependent P0 signal requires explicit wording

The saved results already contain the following estimates for the continuous local nuclear EC cycling module, using the full 21-gene CM score:

| Stage | CM-state adjustment | Standardized coefficient | 95% spatial-block interval |
|---|---|---:|---:|
| P0 | Includes CM cycling and Slit2 | −0.00686 | −0.02236 to +0.00590 |
| P0 | Removes both CM-state covariates | −0.01686 | −0.03294 to −0.00104 |
| P7 | Includes CM cycling and Slit2 | +0.00584 | −0.00419 to +0.01347 |
| P7 | Removes both CM-state covariates | −0.00439 | −0.01370 to +0.00306 |

Removing the spatial/anatomical terms also exposes a small negative P0 coefficient (−0.01587; interval −0.02972 to −0.00090). However, this model still adjusts for CM state. The binary near-versus-other comparison remains uncertain at P0/P7, and the stricter individual-EC CM-marker restriction does not recover a consistent negative full-panel association. The P0 sarcomere-component finding is another exploratory, dimension-specific observation.

The appropriate current statement is: **there is no robust, consistent P0/P7 reduction of the overall maturation panel across the tested definitions; a small P0 relationship depends on the model and deserves further examination.** It would be too strong to say that every model is null or that non-CMs do not influence maturation. This review changes interpretation, not any saved coefficient or cell assignment.

## Priorities for a stronger follow-up

1. Specify the biological question before choosing adjustment: the overall local developmental association versus an association conditional on CM cycling and anatomy. Present a technical-covariate model, then an anatomy-adjusted model, then a CM-state-adjusted model in the same cohort. Existing sensitivity models only partly separate these changes.
2. Validate and report distinct maturation dimensions, alongside the fixed 21-gene panel. A category-balanced or snRNA-validated score should be defined independently of the desired spatial result. Removing one outcome gene does not solve outcome-dependent cohort selection.
3. Compare nearby CMs within reviewed ventricular regions using segmentation-aware adjacency/boundary distances, with sensitivity to nuclear versus whole-cell geometry and tissue barriers. Region matching can complement a smooth XY regression.
4. Harmonize FB and EC analyses on a common CM cohort and explicit count/identity views. Refit FB family/matrix associations with the latest segmentation-supported sources before calling that pathway fully updated. Fit joint measured FB/EC programs where supported; immune/mural states have so far mainly entered these maturation models as broad composition covariates.
5. Quantify reference/segmentation uncertainty and obtain independent spatial sections/animals for biological validation. More bootstrap draws can stabilize numerical intervals but cannot create biological replication or validate an unmeasured sprouting program.

Sources: [maturation and neighborhood construction](../../analysis/map_fibroblast_states_Xenium.py), [reviewed FB refit](../../analysis/refit_fibroblast_Xenium_reviewed.py), [FB regression/prediction](../../analysis/compare_fibroblast_matrix_and_spatial.py), [FB sensitivities](../../analysis/check_fibroblast_spatial_model_sensitivity.py), [EC construction/regression/sensitivities](../../analysis/cardiac_spatial_maturation_EC.py), [all EC coefficients](association_coefficients.tsv), [FB methods](../fibroblast_Xenium_reviewed_refit/association_methods.json), [EC coverage](neighborhood_coverage.tsv), [FB coverage](../fibroblast_Xenium_reviewed_refit/spatial_model_coverage.tsv).
