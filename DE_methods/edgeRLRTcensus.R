# monocle_2.10.1 edgeR_3.24.3
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(monocle))

run_edgeRLRTcensus <- function(data) {

  phenoData <- new("AnnotatedDataFrame",
                   data = data.frame(condition = data$condt,
                                     row.names = colnames(data$tpm)))
  cds <- newCellDataSet(data$tpm,
                        phenoData = phenoData)
  censuscounts <- relative2abs(cds)
  dge <- DGEList(censuscounts, group = data$condt)
  dge <- edgeR::calcNormFactors(dge)
  design <- model.matrix(~data$condt)
  dge <- estimateDisp(dge, design = design)
  fit <- glmFit(dge, design = design)
  lrt <- glmLRT(fit)
  res <- topTags(lrt, n = Inf)

  list(df = data.frame(pval = res$table$PValue,
                       padj = res$table$FDR,
                       row.names = rownames(res$table)))

}