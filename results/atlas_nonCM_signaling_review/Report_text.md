# Relating cardiomyocyte development to signals from surrounding cells

## Main assessment

The heart atlas offers a useful approach that goes beyond comparing which cell types are nearby. It combines spatial proximity, inferred ligand–receptor activity, downstream transcriptional programs and selected perturbation experiments. We should adapt this sequence, while testing whether the molecular activity of neighboring cells adds information beyond anatomy, RNA capture and the cardiomyocyte's own state.

Our latest analysis was narrower: it tested broad and fine neighbor identities and several measured fibroblast programs using conditional regressions and Ridge prediction. Failure of those features to improve the overall 21-gene prediction does not rule out a particular ligand, receptor or nonlinear interaction. It also does not reproduce the atlas's ligand–receptor random-forest model.

A direct performance comparison would additionally require the same CM selection and maturation-score implementation. Our transparent frozen-reference score is not asserted to recover the authors' unavailable precomputed Seurat score, so the two R² values should not be compared as if they evaluated the same task.

## What the paper actually tested

Figure 4e and Extended Data Figure 6 relate the P7 maturation index to graph-based neighborhood levels around multiple non-CM populations. The main figure highlights capillary ECs, cycling ECs and pericytes. The released plotting script selects compact and trabecular ventricular regions, pools the selected CM states, averages maturation within each neighborhood level and smooths those averages. This is an informative description of regional organization. A correlation between averaged levels is not an estimate of a local cellular influence, and graph rings are not fixed distances in micrometers.

Figures 2 and 5 examine communication with COMMOT. The Methods describe a 200-µm range for secreted signals and 50 µm for matrix/contact-associated pairs. COMMOT uses expression and location, including competition among ligand/receptor species, to infer communication scores. Fifty micrometers is a neighborhood allowance, not demonstrated membrane contact. Matrix RNA indicates a potential matrix-producing source; it does not locate the deposited matrix protein.

The paper describes early Postn fibroblasts as communication hubs and capillary ECs as later hubs. Hub status refers to inferred network activity and may depend on source abundance, expression, capture and panel design. It does not directly establish which population promotes maturation. Figure 5e also identifies CMs themselves as major predicted senders to ventricular CMs. Any analysis attributing maturation to non-CMs must therefore keep CM-to-CM signals as a separate comparison.

Figure 5f–h predicts the maturation index from individual-CM COMMOT ligand–receptor features using a random forest with 100 trees and SHAP interpretation. The released script selects COMMOT pair columns and uses a random 80:20 cell split. This can identify nonlinear predictive relationships. Nearby cells can nevertheless occur in both sets, so this is not independent spatial-region or animal validation. High prediction and a large SHAP value do not establish a causal effect or prove that the feature reflects a non-CM source rather than the receiver's expression state.

SigNET and SCENIC nominate downstream CM transcriptional responses. For the highlighted PTPRM example, the paper adds boundary protein staining and gene perturbation: knockout changes cell-cycle and maturation-related expression. This is stronger evidence for a role of the gene than spatial correlation alone. It does not isolate an endothelial source from a CM-intrinsic role or CM-to-CM interactions, and cell-cycle RNA is not a complete measurement of physiological maturation.

## The approach most useful for our biological questions

The main unit should be a defined ventricular CM and its surrounding molecular environment. Analyze each age separately and distinguish working, Ankrd1-associated and other supported CM states. Retain the directly measured 21-gene index, metabolic-regulatory, sarcomere and excitation–contraction components, with cycling separate. A local signal may affect one component without shifting the whole index.

For each sender population, quantify ligand or matrix-program expression around the CM, the relevant receptor expression in the CM, and an independent downstream response. Keep the sender identity attached to the feature: the same ligand produced by a fibroblast, EC, immune cell or neighboring CM should not automatically become one pooled exposure. Use both total local ligand availability and average expression per informative source cell, alongside source counts. This separates more cells from more expression per cell. Keep Col8a1 as one fibroblast family; use continuous Gfpt2/Pi16 and matrix programs alongside uncertain late identities.

The decisive comparison is whether the neighbors add information after the CM's own receptors are known. Compare five nested feature sets with the same learning method: capture/size/CM identity/anatomy/broad context; that baseline plus CM receptors; then CM-to-CM inputs; then sender-resolved non-CM ligands; then ligand–receptor compatibility or interaction features. Run communication inference jointly across cardiac cell types and retain per-edge sender identities, rather than fitting each lineage in isolation. A gain from receptor expression alone supports a receiver-state marker. A reproducible additional gain from local non-CM sender features is stronger evidence of a neighborhood relationship, although still observational.

Use a constrained nonlinear model alongside interpretable regressions. Evaluate all feature sets using the same folds, and keep tuning inside the training data. Hold out contiguous spatial regions rather than randomly selected cells. For an exposure radius of 200 µm, a 150-µm receiver buffer can still allow shared source cells between training and test sets. Where feasible, use a gap of at least twice the maximum source radius and explicitly check source-ID overlap; adapt the validation design to spatial autocorrelation and cell geometry. Keep test-domain communication construction separate when global competition or feature fitting could connect the domains. Spatial validation remains within-section validation until independent hearts are available.

Check supported relationships against matched anatomical regions and region-preserving source-program nulls. Re-run candidate relationships using deposited and core RNA and alternative source confidence definitions. SHAP can describe a well-performing model, but grouped feature ablation and source-specific checks are needed when ligands are correlated. Do not interpret a feature ranking as a list of causal regulators.

## Feasible biological priorities in this dataset

For EC–CM relationships, distinguish vascular proximity from an expressed contact-associated candidate and from the CM's downstream response. The panel includes Ptprm, but not Meis1; therefore, a spatial Meis1 activity map would be inferred rather than directly measured. Oxygen delivery and blood flow remain unmeasured. Spatial RNA cannot choose between perfusion and contact without additional evidence.

For fibroblast–CM relationships, the priority is whether Postn-, Gfpt2/Pi16- or Col8a1-associated surroundings differ in matrix/receptor programs related to CM structural or metabolic development. Fn1, Itga5, Itgb1, Lama4, Lamc1 and Dag1 are measured targets that support a focused feasibility screen. Their presence in the panel does not establish expression, a fully assembled receptor, deposited matrix or signaling activity. At P21, use the measured fibroblast programs across all eligible fibroblasts as well as family labels, because Gfpt2 annotation remains uncertain.

For endocardial–CM relationships, examine the paper's early growth-signal candidates while matching inner-wall anatomy and keeping CM identity explicit. Several relevant genes are measured, but Erbb2 is absent, so some receptor complexes are incompletely observed. For immune and mural sources, begin with supported broad lineages and measured trophic/matrix programs rather than demanding sparse fine labels.

Learn downstream programs from corrected snRNA-seq, test their coherence across libraries, and examine only adequately covered measured targets in Xenium. The published WGCNA modules are unevenly covered: their measured gene counts are 83, 70, 2 and 10. The latter two cannot be treated as complete spatial measurements of the RNA modules. Remove any genes reused in predictor features, retain CM-lineage-supported outcomes, and keep validated projections separate from measured RNA. A NicheNet-style ligand-to-target prioritization could complement SigNET, but would nominate rather than establish signaling mechanisms.

## Interpretation and next analysis

The most informative next analysis is sender-resolved ligand–receptor–response modeling, not another round of niche clustering. The goal is to identify a specific cell source, molecular program, CM response and developmental interval that remain associated after appropriate spatial and receiver-state controls. Existing CellCharter niches can describe anatomy; existing Wormhole matches can compare how validated environments change across stages. Neither should substitute for the local signaling test or be interpreted as cell ancestry.

This review has checked the paper, public scripts and measured-gene feasibility. It has not fitted a new COMMOT/random-forest model or established a new non-CM influence on CM maturation.

## Sources and reproducibility

Wang et al. Spatially guided in vivo single-cell functional genomics of postnatal heart. Nature Cardiovascular Research (2026), DOI 10.1038/s44161-026-00861-z. Relevant sections: Figures 2, 4 and 5; Extended Data Figures 6 and 7; Methods sections 23, 25, 26 and 27. https://pmc.ncbi.nlm.nih.gov/articles/PMC13561863/

Public PIP_seq code was inspected at commit 102d80e6ea5de44f594461c727a4b3c6e1bc3508. Figure4/4.1_Maturation_Model.py lines 101–131 select pair features and make the random split; lines 146–157 fit and evaluate the forest. Figure4/4.2_main_figure_generation.R selects anatomical sources and plots neighborhood-level means. https://github.com/whfbio/PIP_seq/tree/102d80e6ea5de44f594461c727a4b3c6e1bc3508

The released communication examples use a mouse CellChat database, whereas the Methods name CellPhoneDB v4.0. The example SigNET scripts also do not alone recover the full stated Methods workflow. A reproduction should freeze a documented database and input provenance, rather than silently treat these examples as the exact final analysis. The released WebAtlas observation fields do not include the original maturation/COMMOT feature matrices.

Gene and module coverage are audited against the deposited 403-gene target list and Supplementary Table 8. Availability is not detection, specificity or pathway activation. The review is dated 27 October 2026; no original data or earlier report was changed.
