# Readable labeled UMAPs

Larger single-view figures with high-contrast colors, external identity labels, connector lines, subcluster IDs and nucleus counts. Existing coordinates, memberships and annotations are preserved exactly; no geometric separation or UMAP refitting was performed. Nearby or overlapping populations remain nearby or overlapping.

[Download the combined 11-page PDF](all_labeled_umaps.pdf).

Subtype-colored views share a color across clusters with the same proposed identity. Subcluster-colored views give each computational cluster a distinct color. Candidate/unresolved labels and additional review flags remain visible. All identities are provisional RNA-based interpretations. A cluster is not necessarily a separate biological cell type; within-lineage coordinate spaces are independent.

## All cell types

![Reviewed broad identities](all_cell_types_labeled.png)

## Cardiomyocyte

[Subtype view (PNG)](Cardiomyocyte_subtypes_labeled.png) · [Subtype view (PDF)](Cardiomyocyte_subtypes_labeled.pdf) · [Subcluster view (PNG)](Cardiomyocyte_subclusters_labeled.png) · [Subcluster view (PDF)](Cardiomyocyte_subclusters_labeled.pdf)

![Cardiomyocyte: proposed identities](Cardiomyocyte_subtypes_labeled.png)

## Fibroblast

[Subtype view (PNG)](Fibroblast_subtypes_labeled.png) · [Subtype view (PDF)](Fibroblast_subtypes_labeled.pdf) · [Subcluster view (PNG)](Fibroblast_subclusters_labeled.png) · [Subcluster view (PDF)](Fibroblast_subclusters_labeled.pdf)

![Fibroblast: proposed identities](Fibroblast_subtypes_labeled.png)

## Endothelial

[Subtype view (PNG)](Endothelial_subtypes_labeled.png) · [Subtype view (PDF)](Endothelial_subtypes_labeled.pdf) · [Subcluster view (PNG)](Endothelial_subclusters_labeled.png) · [Subcluster view (PDF)](Endothelial_subclusters_labeled.pdf)

![Endothelial: proposed identities](Endothelial_subtypes_labeled.png)

## Mural

[Subtype view (PNG)](Mural_subtypes_labeled.png) · [Subtype view (PDF)](Mural_subtypes_labeled.pdf) · [Subcluster view (PNG)](Mural_subclusters_labeled.png) · [Subcluster view (PDF)](Mural_subclusters_labeled.pdf)

![Mural: proposed identities](Mural_subtypes_labeled.png)

## Immune

[Subtype view (PNG)](Immune_subtypes_labeled.png) · [Subtype view (PDF)](Immune_subtypes_labeled.pdf) · [Subcluster view (PNG)](Immune_subclusters_labeled.png) · [Subcluster view (PDF)](Immune_subclusters_labeled.pdf)

![Immune: proposed identities](Immune_subtypes_labeled.png)

## Reproduce

Run `python analysis/plot_labeled_umaps.py` in the analysis environment. Inputs are the existing annotation CSVs in `results/cell_identities_subtypes`; no expression reanalysis is performed. `manifest.json` lists plotted nucleus counts, coordinate columns and artifact checksums.
