# metagenomeSeq_1.24.1
suppressPackageStartupMessages(library(metagenomeSeq))

run_metagenomeSeq <- function(data) {

  phenoData <- new(
    "AnnotatedDataFrame",
    data = data.frame(
      condition = data$condt,
      row.names = colnames(data$count)))
  obj <- newMRexperiment(
    data$count,
    phenoData = phenoData)
  p <- cumNormStatFast(obj)
  obj <- cumNorm(obj, p = p)
  mod <- model.matrix(~ condition, data = pData(obj))
  res <- fitFeatureModel(obj = obj, mod = mod, coef = 2)
  res <- MRtable(obj = res, number = Inf, by = 2, adjustMethod = "BH")

  list(df = data.frame(pval = res$pvalues,
                       padj = res$adjPvalues,
                       row.names = rownames(res)))

}
