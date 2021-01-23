args = commandArgs(trailingOnly=TRUE)
method <- as.character(args[1]) # the name of evaluation method
work_dir <- as.character(args[2]) # the dir containing all files should not end with "/"
work_dir <- paste0(work_dir, "/")

# set R library
.libPaths(paste0(work_dir, "3.5"))

dataname <- c("GSE74596", "GSE60749-GPL13112", "GSE45719") # all data used for simulation
nsample <- c(50, 100, 500, 1000, 1500, 2000, 2500, 3000, 3500) # number of cells per group
seed <- c(123, 456) # simulation seed
resultpath <- paste0(work_dir, "scDD_simu_result/", method, "/")
savepath <- paste0(work_dir, "scDD_simu_result/comb_results/")

files <- list.files(path = paste0(resultpath), pattern = "rds$")
all_res <- list()
for(nm in dataname) {
  for (sz in nsample) {
    for (sd in seed) {
      filename <- paste0(method, "_", nm, "_", sz, "_", sd, ".rds")
      if (filename %in% files == FALSE) {next()}
      cat(paste0(filename, "\n"))
      res <- readRDS(paste0(resultpath, filename))
      all_res[[paste0(filename)]] <- res
    }
  }
}
saveRDS(all_res, file = paste0(savepath, method, ".rds"))
