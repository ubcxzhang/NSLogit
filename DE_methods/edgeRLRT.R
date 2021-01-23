# edgeR_3.24.3
suppressPackageStartupMessages(library(edgeR))

run_edgeRLRT <- function(data) {

  dge <- DGEList(data$count, group = data$condt)
  dge <- calcNormFactors(dge)
  design <- model.matrix(~data$condt)
  dge <- estimateDisp(dge, design = design)
  fit <- glmFit(dge, design = design)
  lrt <- glmLRT(fit)
  res <- topTags(lrt, n = Inf)

  list(df = data.frame(pval = res$table$PValue,
                       padj = res$table$FDR,
                       row.names = rownames(res$table)))

}