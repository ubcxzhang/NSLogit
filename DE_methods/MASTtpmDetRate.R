# MAST_1.8.2
suppressPackageStartupMessages(library(MAST))

run_MASTtpmDetRate <- function(data) {

  stopifnot(all(names(data$condt) == colnames(data$tpm)))
  grp <- data$condt
  cdr <- scale(colMeans(data$tpm > 0))
  sca <- FromMatrix(exprsArray = log2(data$tpm + 1),
                    cData = data.frame(wellKey = colnames(data$tpm),
                                       grp = grp,
                                       cdr = cdr)) # modified
  zlmdata <- zlm(~cdr + grp, sca) # modified
  res <- lrTest(zlmdata, "grp")

  list(df = data.frame(pval = res[, 'hurdle', 'Pr(>Chisq)'],
                       padj = p.adjust(res[, 'hurdle', 'Pr(>Chisq)'], method="fdr", n=nrow(data$tpm)),
                       row.names = rownames(data$tpm)))

}