# edgeR_3.24.3
suppressPackageStartupMessages(library(edgeR))

run_Wilcoxon <- function(data) {

  tmm <- edgeR::calcNormFactors(data$tpm)
  tpmtmm <- edgeR::cpm(data$tpm, lib.size = tmm * colSums(data$tpm))
  idx <- 1:nrow(tpmtmm)
  names(idx) <- rownames(tpmtmm)
  wilcox_p <- sapply(idx, function(i) {
    wilcox.test(tpmtmm[i, ] ~ data$condt)$p.value
  })

  list(df = data.frame(pval = wilcox_p,
                       padj = p.adjust(wilcox_p, method="fdr", n=nrow(data$tpm)),
                       row.names = names(wilcox_p)))

}