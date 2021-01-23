# DESeq2_1.22.2 monocle_2.10.1 edgeR_3.24.3
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(monocle))

run_monoclecensus <- function(data) {

  phenoData <- new("AnnotatedDataFrame",
                   data = data.frame(condition = data$condt,
                                     row.names = colnames(data$tpm)))
  mon <- newCellDataSet(as.matrix(data$tpm),
                        phenoData = phenoData,
                        expressionFamily = tobit())
  rpc_matrix <- relative2abs(mon)
  mon <- newCellDataSet(cellData = as.matrix(rpc_matrix),
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