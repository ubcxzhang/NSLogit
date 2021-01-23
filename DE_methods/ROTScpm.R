# ROTS_1.10.1 edgeR_3.24.3
suppressPackageStartupMessages(library(ROTS))
suppressPackageStartupMessages(library(edgeR))

run_ROTScpm<- function(data) {
  stopifnot(all(names(data$condt) == colnames(data$count)))
  condition <- as.factor(data$condt)
  levels(condition)[1] <- 1
  levels(condition)[2] <- 2
  names(condition) <- colnames(data$count)
  grp <- as.numeric(condition)

  dge <- DGEList(counts = data$count)
  dge <- edgeR::calcNormFactors(dge)
  cpms <- edgeR::cpm(dge) # modified
  res <- ROTS(data = cpms, groups = grp, B = 1000, K = 1000, log = FALSE, seed = 123)

  list(df = data.frame(
    pval = res$pvalue,
    padj = p.adjust(res$pvalue, method="fdr", n=nrow(data$count)),
    row.names = rownames(res$data)))
}