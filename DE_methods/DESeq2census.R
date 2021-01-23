# DESeq2_1.22.2 monocle_2.10.1
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(monocle))

run_DESeq2census <- function(data) {

  data$condt <- as.character(data$condt)

  phenoData = new("AnnotatedDataFrame",
                  data = data.frame(condition = data$condt,
                                    row.names = colnames(data$tpm)))

  cds <- newCellDataSet(data$tpm,
                        phenoData = phenoData)

  censuscounts <- relative2abs(cds)
  dds <- DESeqDataSetFromMatrix(countData = round(censuscounts),
                                colData = data.frame(condition = data$condt),
                                design = ~condition)
  dds <- DESeq(dds)
  res <- results(dds,
                 contrast = c("condition",
                              levels(factor(data$condt))[1],
                              levels(factor(data$condt))[2]),
                 alpha = 0.05)


  list(df = data.frame(pval = res$pvalue,
                       padj = res$padj,
                       row.names = rownames(res)))

}