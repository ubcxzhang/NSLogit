suppressPackageStartupMessages(library(Seurat))

run_SeuratTobitnofilt <- function(data) {
  tmptpm<- data$tpm
  colnames(tmptpm) <- paste0(data$condt, "__", 1:ncol(data$tpm))
  seur <- CreateSeuratObject(
    tmptpm, project = "scrnaseq",
    names.field = 1, names.delim = "__")
  res <- FindMarkers(seur, ident.1 = levels(factor(data$condt))[1],
    ident.2 = levels(factor(data$condt))[2], test.use = "tobit",
    logfc.threshold = -Inf, min.pct = 0, min.cells = 0)

  df <- data.frame(
    pval = rep(NA, dim(data$count)[1]),
    row.names = rownames(data$count))
  df[rownames(res), "pval"] <- res$p_val
  df[rownames(res), "padj"] <- p.adjust(res$p_val, method = "BH")

  list(df = df)
}