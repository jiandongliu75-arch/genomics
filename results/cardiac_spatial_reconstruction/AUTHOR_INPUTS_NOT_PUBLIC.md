# What prevents exact original Figure 4e replication

Recovered public labels and released code make the analysis rules inspectable, but do not provide all original inputs. The [audit table](author_reconstruction_audit.tsv) distinguishes recovered facts from reconstructed choices.

The released Figure 4 script reads:

- `ref.xenium.{stage}.1011.rds`, containing `Maturation_Features1` and the anatomical `source` field.
- `xenium_p7_metadata_cellneighbour_v2.csv`, containing per-cell neighbor levels.
- External helper/plotting scripts not fully bundled with the figure code.

The public WebAtlas observation fields contain native metrics and `cytoSPACE.final.anno`, but omit `source`, `Maturation_Features1` and the neighbor-level columns. The public browser configuration likewise exposes cell-type sets rather than the original anatomical source assignments. The current author repository tree was checked against commit `102d80e6ea5de44f594461c727a4b3c6e1bc3508`.

For an exact reconstruction, the remaining required artifacts are:

1. The per-cell anatomical source assignments or polygons for the four ventricular compact/trabecular categories, with original barcodes and coordinate conventions.
2. The original maturation-score vector and construction code: assay, normalization, reference population, gene list, control pool, binning, seed and scoring parameters.
3. The per-cell neighbor table or graph plus graph-building code, input cohort, coordinate system, neighborhood parameters and treatment of disconnected/≥10-ring profiles.

Figure 4d uses `v.CM` alone. The Figure 4e code uses multiple CM states within the four anatomical categories. Counts in Figure 4d cannot be used to define Figure 4e membership.

The present release reproduces exact public Figure 4d label counts, the released Figure 4e label rules and R LOESS plotting rule. It reconstructs explicit measured-RNA scores, independent morphology selectors and specified graphs, and recovers the descriptive direction under those definitions. It does **not** assert recovery of the missing original masks, scores or graph, or pixel-exact reproduction of the paper figure. No request has been sent to the authors.
