# limma_3.38.3 edgeR_3.24.3

suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(edgeR))

run_voomlimma <- function(data) {

  dge <- DGEList(data$count, group = data$condt)
  dge <- edgeR::calcNormFactors(dge) # modified
  design <- model.matrix(~data$condt)
  vm <- voom(dge, design = design, plot = FALSE)
  fit <- lmFit(vm, design = design)
  fit <- eBayes(fit)
  res <- topTable(fit, n = Inf, adjust.method = "BH")

  list(df = data.frame(pval = res$P.Value,
                       padj = res$adj.P.Val,
                       row.names = rownames(res)))
}