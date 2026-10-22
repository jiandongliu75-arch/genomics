# Comparison with the paper's released CytoSPACE workflow

**22 October 2026 update:** the public WebAtlas final labels are now accessible. All 403,168 exported cells match the original Xenium barcodes, native centroids and transcript totals. The [published CM-state maps and a new independent default-parameter run](../author_atlas_baseline/README.md) are available separately. Original RNA reference objects, original spatial input-label files and original RNA-to-Xenium pairs remain unavailable; the final exported labels do not verify the granularity of those original inputs. The comparison below describes the preceding independent run.

The preceding independent CytoSPACE analysis is **not an exact reproduction of the authors' pipeline**. Its 259 accepted fine matches are the subset passing additional filters, not the number of cells successfully assigned by CytoSPACE. Native CytoSPACE assigned all 241,449 included segments. Differences below affect the interpretation of the comparison; their individual effects on accuracy have not been isolated experimentally.

The comparison uses the [paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC13561863/) and the authors' [Figure1 scripts](https://github.com/whfbio/PIP_seq/tree/102d80e6ea5de44f594461c727a4b3c6e1bc3508/Figure1), inspected at commit `102d80e6ea5de44f594461c727a4b3c6e1bc3508`.

| Component | Authors' released workflow | Our independent run |
|---|---|---|
| Reference identities | Input exporter selects `celltype.lv2`; their annotation code defines this level using ventricular working, Slit2-associated, proliferating and other CM states, alongside nonmyocytes | All ventricular references share one broad `Ventricular_CM` input label; stage-specific fine groups are attached only after reference matching |
| Xenium identities supplied | Command explicitly includes `--st-cell-type-path` | Broad lineage/chamber candidates are supplied, without ventricular subtype identities |
| Biological resolution | Paper reports seven CM states across the atlas, including atrial/atrioventricular states | Five to nine ventricular groups within each age, 28 stage-specific groups in total |
| Reference expression | Preparation exports the `RNA` count assay; command retains default molecule downsampling | Main run uses DecontX-corrected unnormalized estimates with downsampling off, with an original-count sensitivity run |
| Partition size | 10,000 spots; 16 requested processors | 3,000 spots; four processors, with explicitly randomized input order |
| Reported labels | Postprocessing copies output `CellType` into `cytoSPACE.final.anno` | Native matches are saved separately from filtered family/fine candidates |
| Evaluation | Paper describes manual-annotation comparison and marker-map agreement; released code calculates correlations on common measured genes, including per-cell and per-gene correlations | Adds a 19-gene holdout, repeated seeds, original-count agreement, separate-library agreement and reference-panel cross-library recoverability |

## What is established about the label constraints

The [released run command](https://github.com/whfbio/PIP_seq/blob/102d80e6ea5de44f594461c727a4b3c6e1bc3508/Figure1/1.3_Cytospace_run.sh) **supplies Xenium labels**, even though its output folder is named `singlecell_withoutSTcelltype_*`. The command, rather than that folder name, determines the algorithm's behavior.

The [input-preparation script](https://github.com/whfbio/PIP_seq/blob/102d80e6ea5de44f594461c727a4b3c6e1bc3508/Figure1/1.2_Cytospace.R) sets reference identities to `celltype.lv2`. Their [annotation script](https://github.com/whfbio/PIP_seq/blob/102d80e6ea5de44f594461c727a4b3c6e1bc3508/Figure1/1.0_snRNA_Manual_Annotation.R) defines that level using several CM states. Their [postprocessing script](https://github.com/whfbio/PIP_seq/blob/102d80e6ea5de44f594461c727a4b3c6e1bc3508/Figure1/1.4_Post_Cytospace_Correlation.R) directly uses the output `CellType` as the final spatial annotation.

In upstream CytoSPACE1.1.0 single-cell mode, a supplied spatial label constrains matching to references of that label. Therefore, **if subtype-level spatial labels were supplied, those subtype identities are inputs to transcriptome placement, rather than being newly identified by that placement**. That can be appropriate for reconstructing full reference transcriptomes within previously annotated cell states.

However, the exact contents of the authors' `*_st_cell_type_labels.txt` files have **not been verified here**. The inspected scripts load pre-existing annotated RDS objects and do not provide a complete derivation of these label files. The illustrated reference label levels and the paper's final seven-state naming also need reconciliation against the actual serialized objects. We cannot claim that all final subtype identities were predefined, or assess whether the manual labels used for validation were independent of the input labels, from these scripts alone.

## Why the validation results are not interchangeable

A correlation across hundreds of genes within a matched cell can be positive because the broad lineage profile agrees. This can coexist with low accuracy for predicting the spatial variation of a particular gene across ventricular cells. Averaging by broad cell type further reduces noise. The paper's Extended Data Fig.1d shows cell-type-level correlations, while its released code also includes individual-cell and gene-level correlations.

The released correlation code uses the genes common to measured Xenium and the assigned reference. Those genes were available to the matching algorithm. Our held-out-gene evaluation removes 19 genes from the matching cost and assesses their recovery within ventricular candidates. It tests a different, more demanding prediction question. Requiring fine-group stability across reference libraries and repeated runs is also an additional acceptance rule, not the native CytoSPACE definition of a successful assignment.

Our conservative filter counts should therefore **not be used as an estimate of the authors' mapping yield or evidence that their broad atlas failed**. Conversely, a complete, visually coherent label map does not by itself establish independent recovery of every finer ventricular state.

## Parameter differences and spatial coherence

The authors use 10,000-cell partitions and default reference downsampling; we used 3,000-cell partitions and no molecule downsampling. These are real differences to reproduce. Among our QC-retained reference nuclei, fewer than0.5% per stage exceed the default1,500-molecule cap over the shared panel, so that cap alone is unlikely to explain the difference in accepted fine-group counts. Larger partitions may alter assignments and seed sensitivity; this has not been benchmarked here with the authors' exact labels.

The paper's CytoSPACE methods paragraph describes a cost combining expression similarity and spatial smoothness. However, the released command and the inspected [upstream1.1.0 assignment code](https://github.com/digitalcytometry/cytospace/blob/bb7fe12ab4466ef5e38dcac4c3ece07012060d60/cytospace/cytospace.py) do not show a spatial-neighbor-distance term in the expression assignment cost. Coordinates are attached to matched locations. The paper separately uses **CellCharter** with spatial neighbors for neighborhood analysis. A custom spatial-smoothing modification to CytoSPACE cannot be verified from the inspected sources; it should not be assumed to explain the difference.

## Full transcriptomes and Kcnd2

The authors export full reference RNA matrices and retain the assigned-expression assay. Upstream CytoSPACE nevertheless computes matching on the gene intersection. Supplying a full reference transcriptome permits unmeasured genes to be carried to the selected coordinates; it does not make those genes newly observed by Xenium.

Kcnd2 and Irx5 remain absent from the deposited measured panel. Their reference-derived spatial patterns would be predictions in either workflow. The paper's successful atlas construction does not establish direct validation of an unassayed Kcnd2 gradient.

## Appropriate next comparison

A paper-matched baseline requires the authors' actual reference and Xenium annotation tables or serialized objects, reconciliation of their annotation levels, and their original parameter settings. Its raw maps should first be compared at the authors' CM-state resolution. The additional fine-state stability and gene-holdout checks can then be reported separately, without treating their filtered counts as native alignment yield. No new CytoSPACE fit or claim of exact reproduction is made in this methods comparison.
