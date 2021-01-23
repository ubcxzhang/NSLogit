# edgeR_3.24.3 scran_1.10.2
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(scran))

run_edgeRLRTdeconv <- function(data) {

  condition <- as.factor(data$condt)
  levels(condition)[1] <- 1
  levels(condition)[2] <- 2
  names(condition) <- colnames(data$count)

  sce <- SingleCellExperiment(
    assays = list(counts=data$count),
    colData = data.frame(condition))

  sce <- computeSumFactors(sce, sizes = unique(pmin(c(20, 40, 60, 80, 100), ncol(sce)/2)),
    positive = TRUE)

  dge <- convertTo(sce, type = "edgeR")
  design <- model.matrix(~data$condt)
  dge <- estimateDisp(dge, design = design)
  fit <- glmFit(dge, design = design)
  lrt <- glmLRT(fit)
  tt <- topTags(lrt, n = Inf)

  list(df = data.frame(
    pval = tt$table$PValue,
    padj = tt$table$FDR,
    row.names = rownames(tt$table)))

}