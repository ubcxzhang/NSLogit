
args = commandArgs(trailingOnly=TRUE)
dataname <- args[1]
method <- args[2]
sz <- as.integer(args[3])
i <- as.integer(args[4])
filt <- as.logical(args[5]) # whether to filter data or not
work_dir <- as.character(args[6]) # the dir containing all files should not end with "/"
work_dir <- paste0(work_dir, "/")

# set R library
.libPaths(paste0(work_dir, "3.5"))

cat(paste0(dataname," : ", sz, "/", i, "\n"))
cat(paste0("Filtered: ", filt, "\n"))

method_script <- paste0(work_dir, "MasterProject/DE_methods/")
datapath <- paste0(work_dir, "conquer_data/")
subfilepath <- paste0(work_dir, "conquer_data/")
resultpath <- work_dir

source(paste0(method_script, method, ".R"))
source(paste0(work_dir, "MasterProject/utils.R"))

if (filt == TRUE) {
  savepath <- paste0(resultpath, "sim123_filtered/")
} else {
  savepath <- paste0(resultpath, "sim123/")
}

data <- read_sim_data_from_rds(datapath, dataname)
subsets <- readRDS(paste0(subfilepath, dataname, "_subfile", ".rds"))
# take the subset that used by conquer
subdata <- take_sub(data, subsets, sz, i)
count <- subdata$count; tpm <- subdata$tpm

if (filt == FALSE) {
  count <- count[rowSums(count) > 0, ]
  tpm <- tpm[rownames(count), ]
} else {
  nbr <- 0.25 * ncol(count)
  keep_rows <- rownames(tpm)[which(rowSums(tpm > 1) > nbr)]
  count <- count[match(keep_rows, rownames(count)), ]
  tpm <- tpm[match(keep_rows, rownames(tpm)), ]
}

subdata$count <- count; subdata$tpm <- tpm
res <- get(paste0("run_", method))(subdata)
saveRDS(res, file = paste0(savepath, method, "/", dataname,"_", sz, "_", i, "_", method, ".rds"))

cat("Done \n")

