# DEsingle_1.2.1
suppressPackageStartupMessages(library(DEsingle))

run_DEsingle <- function(data) {

  stopifnot(all(colnames(data$count) == names(data$condt)))
  results <- DEsingle(counts = round(data$count), group = factor(data$condt))
  res <- data.frame(results, stringsAsFactors = FALSE)

  list(df = data.frame(pval = res$pvalue,
                       padj = res$pvalue.adj.FDR,
                       row.names = rownames(res)))

}