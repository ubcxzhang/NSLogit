# samr_3.0 impute_1.56.0
suppressPackageStartupMessages(library(samr))

run_SAMseq <- function(data) {

  SAMseq.test <- SAMseq(round(data$count), as.numeric(as.factor(data$condt)),
    resp.type = "Two class unpaired",
    geneid = rownames(data$count), genenames = rownames(data$count),
    nperms = 100, nresamp = 20, fdr.output = 1)
  SAMseq.result.table <- rbind(SAMseq.test$siggenes.table$genes.up,
    SAMseq.test$siggenes.table$genes.lo)
  SAMseq.FDR <- rep(NA, nrow(data$count))
  SAMseq.FDR[match(SAMseq.result.table[, 1], rownames(data$count))] <-
    as.numeric(SAMseq.result.table[, 5])/100

  list(df = data.frame(
    pval = NA,
    padj = SAMseq.FDR,
    row.names = rownames(data$count)))

}