# Reproduce the released Figure4e plotting rule on explicitly reconstructed
# inputs. This does not claim recovery of the missing author RDS/CSV/masks.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==1)
folder <- args[[1]]
x <- read.delim(file.path(folder,"paper_style_neighbor_level_means.tsv"), check.names=FALSE)
keys <- c("stage","cohort","graph","source","outcome")
groups <- split(seq_len(nrow(x)), interaction(x[keys],drop=TRUE))
result <- list()
for (indices in groups) {
  z <- x[indices,]
  if (nrow(z)<4) next
  # Same default span, degree and family as stats::loess in released R script.
  fit <- tryCatch(suppressWarnings(loess(mean ~ level,data=z)),error=function(e)NULL)
  if (is.null(fit)) next
  grid <- seq(min(z$level),max(z$level),length.out=100)
  pred <- suppressWarnings(predict(fit,newdata=data.frame(level=grid)))
  d <- z[rep(1,length(grid)),keys]
  d$level <- grid;d$LOESS_mean <- pred
  result[[length(result)+1]] <- d
}
write.table(do.call(rbind,result),file.path(folder,"paper_style_R_LOESS_curves.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
writeLines(trimws(capture.output(sessionInfo()),which="right"),file.path(folder,"R_LOESS_session_info.txt"))
