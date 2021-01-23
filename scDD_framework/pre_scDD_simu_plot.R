args = commandArgs(trailingOnly=TRUE)
work_dir <- as.character(args[1]) # the dir containing all files should not end with "/"
work_dir <- paste0(work_dir, "/")

# set R library
.libPaths(paste0(work_dir, "3.5"))

suppressPackageStartupMessages(library(dplyr))

# the dir containing all files
setwd(paste0(work_dir, "scDD_simu_result/comb_results/"))

method <- c("DEsingle", "edgeRLRT", "edgeRLRTcensus",
  "edgeRLRTdeconv", "edgeRLRTrobust", "edgeRQLF", "edgeRQLFDetRate",
  "limmatrend", "MASTcpm", "MASTcpmDetRate", "MASTtpmDetRate", "MASTtpm", "metagenomeSeq",
  "monocle", "monoclecensus", "monoclecount",
  "ROTScpm", "ROTStpm", "ROTSvoom",
  "SAMseq", "scDD", "SeuratBimod", "SeuratBimodIsExpr2",
  "SeuratBimodnofilt", "SeuratTobit", "SeuratTobitnofilt",
  "voomlimma", "Wilcoxon", "ttest",
  "BPSC", "LR", "NsLogit")

nsample <- c(50, 100, 500, 1000, 2500, 1500, 2000, 3000, 3500)
seed <- c(123, 456)
dataname <- c("GSE74596", "GSE60749-GPL13112", "GSE45719")

report <- function(res, gene_type) {
  cur <- data.frame()
  gene_res <- res$df %>% filter(gene==gene_type)
  gene_res$padj[is.na(gene_res$padj)] <- 1 # set padj = NA to 1
  cur[1, "Value"] <- sum(gene_res$padj<0.05)/dim(gene_res)[1] # positive rate
  cur[1, "Time"] <- res$time
  cur[1, "Size"] <- as.character(res$df$size[1])
  cur[1, "Gene"] <- gene_type
  return(cur)
}

report_fdr <- function(res) {
  cur <- data.frame()
  res$df$padj[is.na(res$df$padj)] <- 1
  positive <- res$df %>% filter(padj<0.05)
  FP <- positive %>% filter(gene %in% c("EE", "EP"))
  FDR <- dim(FP)[1]/dim(positive)[1]
  cur[1, "Value"] <- FDR
  cur[1, "Time"] <- res$time
  cur[1, "Size"] <- as.character(res$df$size[1])
  return(cur)
}

report_f1  <- function(res) {
  cur <- data.frame()
  res$df$padj[is.na(res$df$padj)] <- 1
  TP <- res$df %>% filter(padj<0.05&gene %in% c("DE", "DP", "DM", "DB"))
  FN <- res$df %>% filter(padj>=0.05&gene %in% c("DE", "DP", "DM", "DB"))
  FP <- res$df %>% filter(padj<0.05&gene %in% c("EE", "EP"))
  TP <- dim(TP)[1]; FN <- dim(FN)[1]; FP <- dim(FP)[1]
  F1 = 2*TP/(2*TP + FN + FP)
  cur[1, "Value"] <- F1
  cur[1, "Time"] <- res$time
  cur[1, "Size"] <- as.character(res$df$size[1])
  return(cur)
}

report_tpr <- function(res) {
  cur <- data.frame()
  res$df$padj[is.na(res$df$padj)] <- 1
  TP <- res$df %>% filter(padj<0.05&gene %in% c("DE", "DP", "DM", "DB"))
  P  <- res$df %>% filter(gene %in% c("DE", "DP", "DM", "DB"))
  TPR <- dim(TP)[1]/dim(P)[1]
  cur[1, "Value"] <- TPR
  cur[1, "Time"] <- res$time
  cur[1, "Size"] <- as.character(res$df$size[1])
  return(cur)
}

files <- list.files(path = paste0(work_dir, "scDD_simu_result/comb_results/"), pattern = "rds$")
DD <- data.frame(); FDR <- data.frame(); TPR <- data.frame(); F1 <- data.frame()
for (mt in method) {
  mt_file <- paste0(mt, ".rds")
  if (mt_file %in% files == FALSE) {next()}
  res <- readRDS(paste0(mt, ".rds"))
  cat(paste0(mt, "\n"))
  for (exp in names(res)) {
    exp_res <- res[[exp]]
    dd  <- rbind(report(exp_res, "DE"), report(exp_res, "DP"),
      report(exp_res, "DM"), report(exp_res, "DB"))
    fdr <- report_fdr(exp_res)
    tpr <- report_tpr(exp_res)
    f1 <- report_f1(exp_res)
    dd["Method"] <- mt; fdr["Method"] <- mt; tpr["Method"] <- mt
    f1["Method"] <- mt

    DD <- rbind(DD, dd)
    FDR  <- rbind(FDR , fdr); TPR  <- rbind(TPR , tpr);
    F1 <- rbind(F1, f1)
  }
}

scDD_simu_res <- list(FDR=FDR, TPR=TPR, DD=DD, F1=F1)
saveRDS(scDD_simu_res, file = paste0(work_dir, "scDD_simu_result/comb_results/scDD_simu_res.rds"))
cat(paste0("scDD_simu_res.rds saved at", work_dir, "scDD_simu_result/comb_results \n"))
