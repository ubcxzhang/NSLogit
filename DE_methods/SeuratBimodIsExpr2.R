#Seurat_3.1.1
suppressPackageStartupMessages(library(Seurat))

run_SeuratBimodIsExpr2 <- function(data) {
  tmpcount <- data$count
  colnames(tmpcount) <- paste0(data$condt, "__", 1:ncol(data$count))

  seur <- CreateSeuratObject(
    tmpcount, project = "scrnaseq",
    names.field = 1, names.delim = "__",
    is.expr = 2)
  seur <- NormalizeData(seur, normalization.method = "LogNormalize")

  res <- FindMarkers(seur, ident.1 = levels(factor(data$condt))[1],
    ident.2 = levels(factor(data$condt))[2], test.use = "bimod")

  df <- data.frame(
    pval = rep(NA, dim(data$count)[1]),
    row.names = rownames(data$count))
  df[rownames(res), "pval"] <- res$p_val
  df[rownames(res), "padj"] <- p.adjust(res$p_val, method = "BH")

  list(df = df)
}