
get_fdr <- function(res_sub, sz, mt) {
  cur <- data.frame()
  positive <- res_sub %>% filter(padj<0.05)
  FP <- positive %>% filter(status == 0)
  if (nrow(positive) != 0) {
    FDR <- dim(FP)[1]/dim(positive)[1]
  } else {
    FDR <- 0
  }
  cur[1, "Value"] <- FDR
  cur[1, "Size"] <- as.character(sz)
  cur[1, "Method"] <- as.character(mt)
  return(cur)
}

get_tpr <- function(res_sub, sz, mt) {
  cur <- data.frame()
  TP <- res_sub %>% filter(padj<0.05 & status == 1)
  P  <- res_sub %>% filter(status == 1)
  TPR <- dim(TP)[1]/dim(P)[1]
  cur[1, "Value"] <- TPR
  cur[1, "Size"] <- as.character(sz)
  cur[1, "Method"] <- as.character(mt)
  return(cur)
}

get_f1  <- function(res_sub, sz, mt) {
  cur <- data.frame()
  TP <- res_sub %>% filter(padj < 0.05 & status == 1)
  FN <- res_sub %>% filter(padj >= 0.05 & status == 1)
  FP <- res_sub %>% filter(padj < 0.05 & status == 0)
  TP <- dim(TP)[1]; FN <- dim(FN)[1]; FP <- dim(FP)[1]
  F1 = 2*TP/(2*TP + FN + FP)
  cur[1, "Value"] <- F1
  cur[1, "Size"] <- as.character(sz)
  cur[1, "Method"] <- as.character(mt)
  return(cur)
}

get_res <- function(dataname, method, datapath, filtered, path) {
  setwd(datapath)

  if (filtered == TRUE) {file <- "_TPM_1_25p.rds"} else {file <- ".rds"}

  FDR <- data.frame(); TPR <- data.frame(); F1 <- data.frame()

  for (nm in dataname) {
    data <- read_sim_data_from_rds(path, nm)
    subsets <- readRDS(paste0(path, nm, "_subfile", ".rds"))
    truth <- readRDS(paste0(path, nm, "_truth.rds"))

    for (mt in method) {
      resultfile <- paste0(nm, "_", mt, file)
      if (resultfile %in% list.files(path = paste0(datapath), pattern = "rds$") == FALSE) {next()}
      res <- readRDS(paste0(nm, "_", mt, file))

      for (sub in names(res)) {

        if ("df" %in% names(res[[paste0(sub)]]) == FALSE) {next()}

        sz <- strsplit(sub, "\\.")[[1]][2]
        i <- strsplit(sub, "\\.")[[1]][3]
        res_sub <- res[[paste0(sub)]]$df
        cat(paste0(nm, "/", mt, "/", sub, "/", sz, "\n"))

        # take the subset that used by conquer
        subdata <- take_sub(data, subsets, sz, i)
        count <- subdata$count; tpm <- subdata$tpm

        if (filtered == FALSE) {
          count <- count[rowSums(count) > 0, ]
          tpm <- tpm[rownames(count), ]
        } else {
          nbr <- 0.25 * ncol(count)
          keep_rows <- rownames(tpm)[which(rowSums(tpm > 1) > nbr)]
          count <- count[match(keep_rows, rownames(count)), ]
          tpm <- tpm[match(keep_rows, rownames(tpm)), ]
        }

        ngene <- dim(count)[1]
        all_gene <- rownames(count)

        if ("padj" %in% names(res_sub)){ add_padj <- all(is.na(res_sub$padj))} else {add_padj <- TRUE}
        if (add_padj) {res_sub$padj <- p.adjust(res_sub$pval, method="fdr", n=ngene)}

        if ("pval" %in% names(res_sub) == FALSE) {res_sub$pval <- NA}

        res_sub <- res_sub[, c("padj", "pval")]

        omit_gene <- all_gene[all_gene %in% rownames(res_sub) == FALSE]
        omit <- data.frame(
          pval = rep(NA, length(omit_gene)),
          padj = rep(NA, length(omit_gene)),
          row.names = omit_gene)

        res_sub <- rbind(res_sub, omit)
        res_sub$padj[is.na(res_sub$padj)] <- 1
        res_sub$status <- truth[rownames(res_sub), "status"]

        fdr <- get_fdr(res_sub, sz, mt)
        tpr <- get_tpr(res_sub, sz, mt)
        f1 <- get_f1(res_sub, sz, mt)
        FDR <- rbind(FDR, fdr); TPR <- rbind(TPR, tpr); F1 <- rbind(F1, f1)

      }
    }
  }
  return(list(FDR=FDR, TPR=TPR, F1=F1))
}

args = commandArgs(trailingOnly=TRUE)
work_dir <- as.character(args[1]) # the dir containing all files should not end with "/"
work_dir <- paste0(work_dir, "/")

# set R library
.libPaths(paste0(work_dir, "3.5"))

suppressPackageStartupMessages(library(dplyr))

source(paste0(work_dir, "MasterProject/utils.R"))
path <- paste0(work_dir, "conquer_data/") # original data

dataname <- c("GSE74596sim123", "GSE45719sim123", "GSE60749-GPL13112sim123")

method <- c(
  "edgeRLRT",         "edgeRLRTdeconv", "edgeRLRTcensus",  "edgeRQLF",          "edgeRLRTrobust", "edgeRQLFDetRate",
  "ROTSvoom",         "ROTScpm",        "ROTStpm",         "DESeq2",            "DESeq2nofilt",   "DESeq2census",
  "DESeq2betapFALSE", "MASTtpm",        "MASTcpm",         "MASTtpmDetRate",    "MASTcpmDetRate", "monocle",
  "monoclecensus",    "monoclecount",   "SeuratBimod",     "SeuratBimodnofilt", "SeuratTobit",    "SeuratBimodIsExpr2",
  "SAMseq",           "Wilcoxon",       "NODES",           "SCDE",              "BPSC",           "D3E",
  "voomlimma",        "limmatrend",     "metagenomeSeq",   "scDD",              "ttest",          "DEsingle" )

new_method <- c("NsLogit", "LR")

unfilter <- get_res(dataname, new_method, paste0(work_dir, "sim123/"), FALSE, path)
filter <- get_res(dataname, new_method, paste0(work_dir, "sim123_filtered/"), TRUE, path)
conquer_unfilter <- get_res(dataname, method, paste0(work_dir, "conquer_sim123/"), FALSE, path)
conquer_filter <- get_res(dataname, method, paste0(work_dir, "conquer_sim123/"), TRUE, path)

res <- list(
  FDR = rbind(conquer_unfilter$FDR, unfilter$FDR),
  TPR = rbind(conquer_unfilter$TPR, unfilter$TPR),
  F1 = rbind(conquer_unfilter$F1, unfilter$F1))

res_TPM_1_25p <- list(
  FDR = rbind(conquer_filter$FDR, filter$FDR),
  TPR = rbind(conquer_filter$TPR, filter$TPR),
  F1 = rbind(conquer_filter$F1, filter$F1))

saveRDS(res, file = paste0(work_dir, "sim123/sim123_res.rds"))
saveRDS(res_TPM_1_25p, file = paste0(work_dir, "sim123_filtered/sim123_res_TPM_1_25p.rds"))


