# DESeq2_1.22.2 BiocParallel_1.16.6
suppressPackageStartupMessages(library(DESeq2))

run_DESeq2nofilt <- function(data) {

  data$condt <- as.character(data$condt)

  dds <- DESeqDataSetFromMatrix(
    countData = data$count,
    colData = data.frame(condition = data$condt),
    design = ~condition)

  dds <- DESeq(dds, betaPrior = TRUE, minReplicatesForReplace = Inf)

  res <- results(
    dds,
    contrast = c(
      "condition",
      levels(factor(data$condt))[1],
      levels(factor(data$condt))[2]),
    alpha = 0.05,
    cooksCutoff = FALSE,
    independentFiltering = FALSE)

  list(df = data.frame(
    pval = res$pvalue,
    padj = res$padj,
    row.names = rownames(res)))
}