# edgeR_3.24.3
suppressPackageStartupMessages(library(edgeR))

run_edgeRQLFDetRate <- function(data) {

  dge <- DGEList(data$count, group = data$condt)
  dge <- calcNormFactors(dge)
  cdr <- scale(colMeans(data$count > 0))
  design <- model.matrix(~cdr + data$condt)
  dge <- estimateDisp(dge, design = design)
  fit <- glmQLFit(dge, design = design)
  lrt <- glmQLFTest(fit)
  res <- topTags(lrt, n = Inf)


  list(df = data.frame(pval = res$table$PValue,
                       padj = res$table$FDR,
                       row.names = rownames(res$table)))

}