args = commandArgs(trailingOnly=TRUE)
work_dir <- as.character(args[1]) # the dir containing all files should not end with "/"
work_dir <- paste0(work_dir, "/")

method <- c(
  "BPSC", "DEsingle", "edgeRLRT", "edgeRLRTdeconv",
  "edgeRLRTrobust", "edgeRQLF", "edgeRQLFDetRate", "limmatrend",
  "LR", "MASTcpm", "MASTcpmDetRate", "metagenomeSeq",
  "monoclecount", "NsLogit", "ROTScpm", "ROTSvoom",
  "SAMseq", "scDD", "SeuratBimod", "SeuratBimodIsExpr2",
  "SeuratBimodnofilt", "voomlimma")

setwd(paste0(work_dir, "GSE135893_file/"))

# set R library
.libPaths(paste0(work_dir, "3.5"))

suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(pheatmap))

paper_file <- paste0(work_dir, "GSE135893_file/Disease_vs_Control/")
se <- read.csv(file = paste0(paper_file, "Differentiating_Ciliated_disease_vs_control_.csv"),  header = TRUE)
se <- as.data.frame(se)
rownames(se) <- se$X

delist <- list()
for (mt in method) {
  res <- readRDS(paste0(mt, "_count_result/Differentiating_Ciliated.rds"))
  if (all(rownames(se) %in% rownames(res$df))) {
    res <- res$df[rownames(se), ]
  } else {
    common_gene <- rownames(se)[rownames(se) %in% rownames(res$df)]
    res <- res$df[common_gene, ]
  }
  cat(paste0(mt, ": Finish\n"))
  res[is.na(res$padj), "padj"] <- 1
  delist[[paste0(mt)]] <- rownames(res[res$padj < 0.05, ])
}
se$padj <- se$p_val_adj
se[is.na(se$padj), "padj"] <- 1
delist[["paper"]] <- rownames(se[se$padj < 0.05, ])

common <- matrix(0, nrow = length(delist), ncol = length(delist))
rownames(common) <- names(delist)
colnames(common) <- names(delist)

for (mt in names(delist)) {
  common[paste0(mt), paste0(mt)] <- length(delist[[paste0(mt)]])
  res1 <- delist[[paste0(mt)]]
  for (mt2 in names(delist)) {
    if (mt2 == mt) {next()}
    res2 <- delist[[paste0(mt2)]]
    common[paste0(mt), paste0(mt2)] <- as.integer(sum(res2 %in% res1))
  }
}

p <- pheatmap(common,
  display_numbers = T, number_color = "black", number_format = "%i",
  color = colorRampPalette(brewer.pal(11, "Spectral"))(30),
  cluster_rows = F, cluster_cols = F, fontsize_number = 5)

png(filename = paste0(work_dir,"Figure/Differentiating_Ciliated.png"), width = 5.5, height = 5, units = 'in', res = 300)
p
dev.off()

cat(paste0("Differentiating_Ciliated.png is saved at", work_dir,"Figure \n"))
