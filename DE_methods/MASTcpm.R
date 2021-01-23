# MAST_1.8.2 edgeR_3.24.3
suppressPackageStartupMessages(library(MAST))
suppressPackageStartupMessages(library(edgeR))

run_MASTcpm <- function(data) {

  stopifnot(all(names(data$condt) == colnames(data$count)))
  grp <- data$condt
  names(grp) <- colnames(data$condt)
  dge <- DGEList(counts = data$count)
  dge <- edgeR::calcNormFactors(dge)
  cpms <- edgeR::cpm(dge) # modified
  sca <- FromMatrix(exprsArray = log2(cpms + 1),
    cData = data.frame(wellKey = colnames(data$count),
      grp = grp))
  zlmdata <- zlm(~grp, sca) # modified
  res <- lrTest(zlmdata, "grp")

  pval =  res[, 'hurdle', 'Pr(>Chisq)']
  list(df = data.frame(pval = pval,
      padj = p.adjust(pval, method="fdr", n=nrow(data$count)),
      row.names = rownames(data$count)))
}