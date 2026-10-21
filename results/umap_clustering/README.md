# snRNA-seq UMAP clustering: GSE290989

73,052 QC-retained mouse heart nuclei from P0, P7, P14 and P21 (two libraries per age).

![Leiden clusters and ages](umap_leiden_clusters.png)

## Method

Leiden community detection was run on the existing expression-space 15-neighbor graph, using the igraph implementation, seed 42 and iterations until convergence. The working resolution is 0.8. Resolutions 0.4 and 1.2 provide a sensitivity comparison; the resolution is not claimed to be optimal.

The underlying graph uses 30 PCs of scaled, log-normalized DecontX-corrected expression, selected from 2,500 variable genes across libraries with mitochondrial genes excluded. UMAP was computed previously (min_dist 0.3, 200 epochs, seed 42). Clustering uses the expression-neighbor graph rather than distances in the two-dimensional UMAP. No batch integration or age regression was applied. Age-associated differences can therefore contribute to clustering.

Cluster IDs are computational labels, not cell-type names. Existing preliminary broad lineage annotations are retained separately. `cluster_summary.tsv` flags clusters with fewer than 50 nuclei or less than 80% agreement with a single preliminary lineage. Such clusters require review; neither small size nor mixed labels establish a novel population.

## Published data

| File | Contents |
|---|---|
| `umap_leiden_clusters.png` / `.pdf` | Primary clusters with numeric labels, alongside the same nuclei colored by age |
| `clustering_resolution_comparison.png` / `.pdf` | Resolutions 0.4, 0.8 and 1.2 |
| `cell_clusters_and_umap.csv.gz` | Every retained nucleus: unique ID, barcode, sample, age, replicate, UMAP1/2, cluster IDs at all resolutions, preliminary lineage, and original/corrected QC metrics |
| `cluster_summary.tsv` | Cluster size, dominant preliminary lineage, lineage agreement, QC summaries, descriptive markers and review flags |
| `cluster_counts_by_age.tsv` | Nucleus counts by cluster and age |
| `cluster_counts_by_sample.tsv` | Nucleus counts by cluster and library |
| `cluster_counts_by_preliminary_lineage.tsv` | Cross-tabulation against preliminary broad annotations |
| `cluster_top_markers.tsv` | Up to 30 descriptive markers per cluster; not a significance-tested differential-expression result |
| `cluster_pseudobulk_corrected_counts.tsv.gz` | All 32,285 input genes summed within each cluster; fractional corrected counts, not rounded or normalized |
| `lineage_marker_corrected_counts.csv.gz` | Per-nucleus corrected counts for the canonical broad-lineage marker panel |
| `umap_clusters_metadata_only.h5ad` | Scanpy-readable nucleus metadata, cluster assignments, PCA/UMAP coordinates and neighbor graph; **no expression matrix** |
| `requirements.txt` | Pinned top-level Python dependencies |
| `manifest.json` | Input SHA-256, parameters, software versions and published-file checksums |

The complete cell-by-gene corrected/raw count matrix is approximately 1.2 GB and is not included in ordinary Git files. A full annotated copy is prepared in the analysis workspace as `snrna_qc_clustered.h5ad`. GitHub Release upload currently requires network access to `uploads.github.com`; these repository files do not imply the full expression matrix has been uploaded.

Marker ranking uses the difference in mean log1p(normalized counts per 10,000) between a cluster and other nuclei, requiring detection in at least 20% of the cluster and a detection-fraction difference of at least 10 percentage points. Mitochondrial and ribosomal genes are excluded from ranked candidates. Rankings are descriptive and are not replicate-aware differential-expression tests. Corrected fractional nonzero counts count as detected.

## Read the published data

```python
import pandas as pd
import scanpy as sc

cells = pd.read_csv('cell_clusters_and_umap.csv.gz', index_col=0)
clusters = pd.read_csv('cluster_summary.tsv', sep='\t')
metadata = sc.read_h5ad('umap_clusters_metadata_only.h5ad')
# metadata.X has zero genes: this file is for coordinates/metadata, not expression analysis.
counts = pd.read_csv('cluster_pseudobulk_corrected_counts.tsv.gz', sep='\t', index_col=0)
```

Reproduce from the original full QC-filtered H5AD and the published annotation table:

```bash
python -m pip install -r results/umap_clustering/requirements.txt
python analysis/umap_clustering.py \
  --input /path/to/snrna_qc_filtered.h5ad \
  --annotations results/umap_clustering/cell_clusters_and_umap.csv.gz \
  --outdir /path/to/new_output
```

Original source matrices: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE290989

The full input H5AD is independently QC-filtered, doublet-screened with scDblFinder and ambient-corrected with DecontX. The original source matrices alone do not reproduce that input without the preceding QC workflow. Empty-droplet background was unavailable. Nucleus sampling and differential QC losses (especially at P7) affect composition comparisons. Earlier reports in the repository describe QC and preliminary broad-lineage profiles.
