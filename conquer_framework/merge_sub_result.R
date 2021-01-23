
# Run this script to merge result file for each size and replicates
args = commandArgs(trailingOnly=TRUE)
dataname <- as.character(args[1])
method <- as.character(args[2])
filt <- as.logical(args[3])
work_dir <- as.character(args[4]) # the dir containing all files should not end with "/"
work_dir <- paste0(work_dir, "/")
resultpath <- work_dir

# set R library
.libPaths(paste0(work_dir, "3.5"))

if (filt == TRUE) {
  savepath <- savepath <- paste0(resultpath, "sim123_filtered/")
} else {
  savepath <- paste0(resultpath, "sim123/")
}

setwd(paste0(savepath, method))

rep <- c(1,5,5,5)
if (dataname == "GSE60749-GPL13112sim123") {
  sizes <- c(90, 48, 24, 12)
} else if (dataname == "GSE45719sim123") {
  sizes <- c(50, 24, 12, 6)
} else if (dataname == "GSE74596sim123") {
  sizes <- c(44, 22, 12, 6)
}

files <- list.files(pattern = "rds$")
res <- list(); j <- 1
for (sz in sizes) {
  for (i in 1:rep[j]) {
    if (paste0(dataname, "_", sz, "_", i, "_", method, ".rds") %in% files == FALSE) {next()}
    res[[paste0(method,'.',sz,'.',i)]] <- readRDS(paste0(dataname, "_", sz, "_", i, "_", method, ".rds"))
  }
  j <- j + 1
}

if (filt == TRUE) {
  saveRDS(res, file=paste0(savepath, dataname, "_", method, "_TPM_1_25p.rds"))
} else {
  saveRDS(res, file=paste0(savepath, dataname, "_", method, ".rds"))
}

# Rscript merge_sub_result.R "GSE45719sim123" "NsLogit" "FALSE"
# Rscript merge_sub_result.R "GSE60749-GPL13112sim123" "NsLogit" "FALSE"
# Rscript merge_sub_result.R "GSE74596sim123" "NsLogit" "FALSE"
# sbatch experiment.sh "GSE74596sim123" "LR" 12 4 "TRUE"
# sbatch experiment.sh "GSE74596sim123" "LR" 12 4 "FALSE"
