# scDD_1.6.1, original conquer script is not working
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(scDD))
suppressPackageStartupMessages(library(scran))

run_scDD <- function(data) {

  # condition should be 1 or 2
  condition <- as.factor(data$condt)
  levels(condition)[1] <- 1
  levels(condition)[2] <- 2
  names(condition) <- colnames(data$count)

  sce <- SingleCellExperiment(
    assays = list(counts=as.matrix(data$count)), # the normcounts assays slot contains one matrix
    colData = data.frame(condition))

  datNorm.scran <- scDD::preprocess(sce,
                                    zero.thresh = 1,
                                    median_norm = TRUE)

  condition <- condition[colnames(datNorm.scran)]
  names(condition) <- colnames(datNorm.scran)

  SDSumExp <- SingleCellExperiment(
    assays = list(normcounts = normcounts(datNorm.scran)),
    colData = data.frame(condition))

  prior_param <- list(alpha = 0.01, mu0 = 0, s0 = 0.01, a0 = 0.01, b0 = 0.01)

  scd <- scDD(SDSumExp, prior_param = prior_param, testZeroes = FALSE,
    param = BiocParallel::MulticoreParam(workers = 1),
    condition = "condition", min.size = 3, min.nonzero = NULL,
    categorize = FALSE)

  res <- results(scd)

  list(df = data.frame(
    pval = res$nonzero.pvalue,
    padj = res$nonzero.pvalue.adj,
    row.names = rownames(res)))
}
