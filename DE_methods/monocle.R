# monocle_2.10.1
suppressPackageStartupMessages(library(monocle))

run_monocle <- function(data) {

  phenoData <- new("AnnotatedDataFrame",
                   data = data.frame(condition = data$condt,
                                     row.names = colnames(data$tpm)))
  mon <- newCellDataSet(as.matrix(data$tpm),
                        phenoData = phenoData,
                        expressionFamily = tobit())
  res <- differentialGeneTest(mon, fullModelFormulaStr = "~condition")


  list(df = data.frame(pval = res$pval,
                       padj = res$qval,
                       row.names = rownames(res)))

}