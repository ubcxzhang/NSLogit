# ROTS_1.10.1
suppressPackageStartupMessages(library(ROTS))

run_ROTStpm <- function(data) {

  stopifnot(all(names(data$condt) == colnames(data$tpm)))
  condition <- as.factor(data$condt)
  levels(condition)[1] <- 1
  levels(condition)[2] <- 2
  names(condition) <- colnames(data$count)
  grp <- as.numeric(condition)
  res <- ROTS(data = data$tpm, groups = grp, B = 1000, K = 1000, log = FALSE, seed = 123)

  list(df = data.frame(pval = res$pvalue,
                       padj = p.adjust(res$pvalue, method="fdr", n=nrow(data$tpm)),
                       row.names = rownames(data$tpm)))
}