# Segmentation in the Heart Atlas: method, evidence and implications

## Assessment

The authors chose a scientifically appropriate segmentation strategy and performed a relevant comparison with expert-drawn boundaries. Their evidence supports improvement over nucleus expansion. However, absolute boundary agreement remains modest and substantially lower at P21 in the published benchmark. The study does not establish uniformly reliable cell boundaries or transcript ownership at the interfaces most relevant to our subtype analysis.

## What the authors did

The Xenium assay used 10-µm frozen heart sections and the Xenium Cell Segmentation Staining Kit. DAPI marked nuclei; ATP1A1, CD45 and E-cadherin marked boundaries; 18S rRNA marked cell interiors; and α-SMA/vimentin provided interior protein information. The on-instrument pipeline generated cell matrices, transcript coordinates/quality scores and boundary polygons. Thus, their main segmentation was not simply uniform expansion around each nucleus.

The Methods section “Quality control and benchmarking of spatial segmentation” describes 10 randomly selected fields spanning early and late hearts and compact versus trabecular/vascular-rich contexts. Expert curators manually outlined cells. Automated masks were matched to manual masks at an overlap threshold of at least 90%. Unmatched automated masks counted as false positives; unmatched manual masks counted as false negatives. A dedicated evaluation considered multinucleated CMs.

Supplementary Figure 1b–e reports the comparison at P0 and P21, including profiles with no visible nucleus and profiles with more than one nucleus. The legend specifies five images per age, except the P21 nucleus-expansion multi-nucleus category, which has three images. The reported n counts image fields; independent-heart replication of this benchmark is not established by the legend. A nucleus-free CM profile in a thin section can be a valid cross-section of a nucleated cell; it should not automatically be removed.

## What the benchmark supports

Membrane and interior staining provide information that nucleus expansion lacks, particularly for large, elongated or multinucleated CMs and for sectioned cytoplasm without a visible nucleus. The manual comparison explicitly targets an important cardiac failure mode. Supplementary Figure 1 shows multimodal segmentation outperforming nucleus expansion in the displayed comparisons.

The absolute scores qualify that conclusion. Visual readings from panel b place the multimodal score at approximately 45–50% at P0 and around 10% at P21. These are approximate plot readings, not recovered source-data values. The stringent 90% overlap requirement makes this a demanding boundary test. Nevertheless, the marked late-stage deterioration means that relative improvement cannot be presented as consistently accurate segmentation throughout postnatal development.

## Limits of the evidence

The reported score is TP/(TP+FP+FN). The Methods call it an accuracy/F1-style metric and the supplement calls it average precision. It is an object-level Jaccard-style detection score at a single overlap threshold; it is not the standard F1 score, 2TP/(2TP+FP+FN), or a precision–recall-area average precision. Such single-threshold scores are used in segmentation benchmarks, and the naming does not invalidate the comparison. The formula remains interpretable, but the terminology and exact overlap definition should be clarified for reproducibility.

The text does not specify whether “90% overlap” is intersection-over-union or another normalization, and does not fully report matching rules, curator blinding, inter-curator agreement or field-level masks sufficient to reconstruct the benchmark. The on-instrument segmentation software version and relevant parameters are also not given in the cited sample-preparation/benchmark sections. The benchmark is useful but not a complete independent anatomical ground truth.

Boundary overlap and transcript ownership are different outcomes. The figure does not quantify transcript reassignment specifically at CM–FB, CM–EC or peri-arterial interfaces, nor the resulting stability of fine labels and maturation–neighbor associations. A low object score does not mean the same fraction of transcripts are wrong: a mask can fail a strict overlap threshold while retaining many correctly assigned interior molecules. Conversely, a small geometric boundary error can matter greatly for a thin EC or fibroblast beside a transcript-rich CM.

Agreement with snRNA-seq, CytoSPACE placement, cell-type correlations or selected marker images cannot by themselves resolve this issue. Broad expression agreement can coexist with systematic local mixing, and matched/imputed expression is not independent transcript-boundary validation. Neither the publication nor our conservative resegmentation establishes that one particular computational assignment is universally correct.

## Consequences for our analysis

Retain the authors' multimodal boundaries as a biologically informed baseline. Use whole-cell evidence for subtype identification where supported; add nuclear evidence as corroboration and as a check for clear conflicts, rather than requiring detection of every subtype marker inside the nucleus. Retain valid nucleus-free CM profiles when membrane/interior evidence supports them. Nuclear absence and a supported alternative nuclear identity are distinct situations.

Repeat our key findings using the original boundaries, a boundary-core transcript layer and nuclear-supported sensitivities. Compare identical focal cells where possible. Preserve uncertainty at thin-cell and peri-arterial interfaces, and explicitly examine whether P21 label loss tracks geometry, capture or marker/reference disagreement. P21 Gfpt2 disagreement already occurs before nuclear confirmation, so segmentation cannot currently be identified as its sole cause.

For a stronger validation, manually adjudicate sampled CM–FB, CM–EC and vessel-wall interfaces across ages, with explicit cell-type-aware split/merge errors and transcript reassignment rates. Reassess subtype proportions and CM outcome associations under controlled boundary erosion/dilation and alternative segmentations. Nuclear-supported molecules are valuable anchors, but cannot supply the complete transcriptional or geometric definition of a large CM.

## Sources

Wang H et al. Spatially guided in vivo single-cell functional genomics of postnatal heart. Nature Cardiovascular Research (2026). Methods: “Xenium sample preparation and Xenium instrument run” and “Quality control and benchmarking of spatial segmentation”; Figure 1b; Supplementary Figure 1b–e. https://doi.org/10.1038/s44161-026-00861-z

Public article: https://pmc.ncbi.nlm.nih.gov/articles/PMC13561863/

Supplementary information inspected directly: https://static-content.springer-cdn.com/esm/art%3A10.1038%2Fs44161-026-00861-z/MediaObjects/44161_2026_861_MOESM1_ESM.pdf . Retrieved through the Europe PMC supplementary-files service; PDF pages 3–4 contain Supplementary Figure 1 and its legend. Journal PDF not redistributed.
