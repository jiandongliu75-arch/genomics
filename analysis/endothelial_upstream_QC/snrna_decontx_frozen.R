args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 1L) stop('Usage: Rscript snrna_decontx.R SAMPLE_DIRECTORY')
out <- normalizePath(args[1])
suppressPackageStartupMessages({
    library(Matrix)
    library(SingleCellExperiment)
    library(decontX)
    library(scDblFinder)
})
set.seed(42)
x <- as(readMM(gzfile(file.path(out, 'decontx_input.mtx.gz'))), 'dgCMatrix')
cells <- read.delim(file.path(out, 'decontx_cells.tsv'), stringsAsFactors=FALSE)
stopifnot(ncol(x) == nrow(cells), all(x@x >= 0))
colnames(x) <- cells$cell_id
sce <- SingleCellExperiment(assays=list(counts=x))
rownames(sce) <- sprintf('gene_%d', seq_len(nrow(x)))
z <- as.integer(factor(cells$cluster))
# scDblFinder is primary; Scrublet automatic calls remain diagnostic only.
sce <- scDblFinder(sce, clusters=z, verbose=TRUE)
doublet_metrics <- data.frame(cell_id=colnames(sce),
    scdblfinder_score=sce$scDblFinder.score, scdblfinder_class=sce$scDblFinder.class)
write.csv(doublet_metrics, file.path(out, 'scdblfinder_metrics.csv'), row.names=FALSE)
keep <- sce$scDblFinder.class == 'singlet'
sce <- sce[, keep]
z <- z[keep]
x <- counts(sce)
fit_args <- list(x=sce, z=z, maxIter=500, seed=12345L, verbose=TRUE)
sce <- do.call(decontX, fit_args)
corrected <- assay(sce, 'decontXcounts')
stopifnot(identical(dim(corrected), dim(x)), all(corrected@x >= 0))
# Sparse interchange avoids the expensive decimal formatting of millions of entries.
# Slots are documented little-endian CSC arrays; final user files are standard H5AD.
writeBin(corrected@x, file.path(out, 'decontx_values.f64'), size=8, endian='little')
writeBin(corrected@i, file.path(out, 'decontx_indices.i32'), size=4, endian='little')
writeBin(corrected@p, file.path(out, 'decontx_indptr.i32'), size=4, endian='little')
write.table(matrix(dim(corrected), nrow=1), file.path(out, 'decontx_shape.txt'),
    row.names=FALSE, col.names=FALSE)
metrics <- data.frame(cell_id=colnames(sce),
    decontx_contamination=sce$decontX_contamination,
    raw_umis=Matrix::colSums(x), corrected_umis=Matrix::colSums(corrected))
write.csv(metrics, file.path(out, 'decontx_metrics.csv'), row.names=FALSE)
capture.output(str(metadata(sce)), file=file.path(out, 'decontx_model_metadata.txt'))
saveRDS(metadata(sce), file.path(out, 'decontx_model_metadata.rds'))
capture.output(sessionInfo(), file=file.path(out, 'decontx_R_sessionInfo.txt'))
cat('DECONTX_DONE', basename(out), 'nuclei', ncol(x),
    'median_contamination', median(sce$decontX_contamination), '\n')
