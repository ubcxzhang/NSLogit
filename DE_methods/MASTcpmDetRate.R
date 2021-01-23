# MAST_1.8.2 edgeR_3.24.3
suppressPackageStartupMessages(library(MAST))
suppressPackageStartupMessages(library(edgeR))

run_MASTcpmDetRate <- function(data) {

  stopifnot(all(names(data$condt) == colnames(data$count)))
  grp <- data$condt
  cdr <- scale(colMeans(data$count > 0))
  dge <- DGEList(counts = data$count)
  dge <- edgeR::calcNormFactors(dge)
  cpms <- edgeR::cpm(dge) # modified
  sca <- FromMatrix(exprsArray = log2(cpms + 1),
                    cData = data.frame(wellKey = colnames(data$count),
                                       grp = grp,
                                       cdr = cdr))# modified
  zlmdata <- zlm(~cdr + grp, sca) # modified
  res <- lrTest(zlmdata, "grp")

  list(df = data.frame(pval = res[, 'hurdle', 'Pr(>Chisq)'],
                       padj = p.adjust(res[, 'hurdle', 'Pr(>Chisq)'], method="fdr", n=nrow(data$count)),
                       row.names = rownames(data$count)))

}