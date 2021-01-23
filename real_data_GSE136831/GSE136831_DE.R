args = commandArgs(trailingOnly=TRUE)
filename <- as.character(args[1])
method <- as.character(args[2])
work_dir <- as.character(args[3]) # the dir containing all files should not end with "/"
work_dir <- paste0(work_dir, "/")

# set R library
.libPaths(paste0(work_dir, "3.5"))
library(Matrix)

method_script <- paste0(work_dir, "MasterProject/DE_methods/")
realdata_part_file <- paste0(work_dir, "GSE136831_count/")
realdata_file <- paste0(work_dir, "GSE136831_file/")

source(paste0(method_script, method, ".R"))
sct <- readRDS(paste0(realdata_part_file, filename, ".rds"))
sct$count <- as.matrix(sct$count)

# remove COPD
sct$condt <- as.character(sct$condt)
names(sct$condt) <- colnames(sct$count)
select_cell <- names(sct$condt[sct$condt != "COPD"])
sct$condt <- sct$condt[select_cell]
sct$count <- sct$count[ ,select_cell]

# apply de
time <- system.time({res <- get(paste0("run_", method))(sct)})
res$time <- time[3]
res$session_info <- sessionInfo() # session info
saveRDS(res, file = paste0(realdata_file, method,"_count_result/", filename, ".rds"))
