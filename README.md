# Postnatal heart snRNA-seq analysis

Independent QC and exploratory clustering of the mouse-heart dataset [GSE290989](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE290989), covering P0, P7, P14 and P21 with two libraries per age.

- [UMAP figures, cluster assignments and downloadable analysis tables](results/umap_clustering/README.md)
- [QC report (PDF)](snRNA_QC_report.pdf)
- [Preliminary lineage profiles (PDF)](heart_lineage_profiles.pdf)
- [Reproducible clustering script](analysis/umap_clustering.py)

![Leiden clusters and ages](results/umap_clustering/umap_leiden_clusters.png)

73,052 of 98,371 input nuclei passed independent QC, doublet screening and ambient correction. The primary Leiden resolution of 0.8 yields 24 computational clusters; resolutions 0.4 and 1.2 yield 15 and 29 clusters. The UMAP is not batch-integrated, and cluster IDs are not final biological subtype annotations.

The repository includes per-nucleus UMAP coordinates and metadata, a metadata-only H5AD, lineage-marker expression, and cluster-level pseudobulk counts for all 32,285 input genes. The full 1.2 GB cell-by-gene expression H5AD is prepared in the analysis workspace but has not been uploaded: GitHub Release asset uploads require access to `uploads.github.com`.
