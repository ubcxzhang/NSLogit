# limma_3.38.3 edgeR_3.24.3
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(edgeR))

run_limmatrend <- function(data) {
  dge <- DGEList(data$count, group = data$condt)
  dge <- edgeR::calcNormFactors(dge) # modified
  design <- model.matrix(~data$condt)
  y <- new("EList")
  y$E <- edgeR::cpm(dge, log = TRUE, prior.count = 3)
  fit <- lmFit(y, design = design)
  fit <- eBayes(fit, trend = TRUE, robust = TRUE)
  res <- topTable(fit, n = Inf, adjust.method = "BH")

  list(df = data.frame(pval = res$P.Value,
      padj = res$adj.P.Val,
      row.names = rownames(res)))
}