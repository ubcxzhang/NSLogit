# monocle_2.10.1
suppressPackageStartupMessages(library(monocle))

run_monoclecount <- function(data) {

  phenoData <- new(
    "AnnotatedDataFrame",
    data = data.frame(
      condition = data$condt,
      row.names = colnames(data$count)))

  mon <- newCellDataSet(
    as.matrix(data$count),
    phenoData = phenoData,
    lowerDetectionLimit = 0.5,
    expressionFamily = negbinomial.size())

  mon <- estimateSizeFactors(mon)
  mon <- estimateDispersions(mon)
  res <- differentialGeneTest(mon, fullModelFormulaStr = " ~ condition")

  list(df = data.frame(pval = res$pval,
                       padj = res$qval,
                       row.names = rownames(res)))

}