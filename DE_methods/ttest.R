# edgeR_3.24.3
suppressPackageStartupMessages(library(edgeR))

run_ttest <- function(data) {

  tmm <- edgeR::calcNormFactors(data$tpm)
  tpmtmm <- edgeR::cpm(data$tpm, lib.size = tmm * colSums(data$tpm))
  logtpm <- log2(tpmtmm + 1)
  idx <- seq_len(nrow(logtpm))
  names(idx) <- rownames(logtpm)
  ttest_p <- sapply(idx, function(i) {
    t.test(logtpm[i, ] ~ data$condt)$p.value
  })

  list(df = data.frame(pval = ttest_p,
                       padj = p.adjust(ttest_p, method="fdr", n=nrow(data$tpm)),
                       row.names = rownames(data$tpm)))

}