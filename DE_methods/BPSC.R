suppressPackageStartupMessages(library(BPSC))
suppressPackageStartupMessages(library(edgeR))

run_BPSC <- function(data) {
  cpms <- edgeR::cpm(data$count, lib.size = colSums(data$count) * edgeR::calcNormFactors(data$count))
    controlIds <- which(data$condt == levels(factor(data$condt))[1])
    design <- model.matrix(~ data$condt)
    coef <- 2
    resbp <- BPglm(data = cpms, controlIds = controlIds,
      design = design, coef = coef, estIntPar = FALSE)

  list(df = data.frame(
    pval = resbp$PVAL,
    padj = p.adjust(resbp$PVAL, method="fdr", n=nrow(data$count)),
    row.names = names(resbp$PVAL)))
}