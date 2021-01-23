# MAST_1.8.2
suppressPackageStartupMessages(library(MAST))

run_MASTtpm <- function(data) {

  stopifnot(all(names(data$condt) == colnames(data$tpm)))
  grp <- data$condt
  sca <- FromMatrix(exprsArray = log2(data$tpm + 1),
                    cData = data.frame(wellKey = colnames(data$tpm),
                                       grp = grp)) # modified
  zlmdata <- zlm(~grp, sca) # modified
  res <- lrTest(zlmdata, "grp")

  list(df = data.frame(pval = res[, 'hurdle', 'Pr(>Chisq)'],
                       padj = p.adjust(res[, 'hurdle', 'Pr(>Chisq)'], method="fdr", n=nrow(data$tpm)),
                       row.names = rownames(data$tpm)))

}