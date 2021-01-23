
# Run DE on nm_nsample_seed.rds
args = commandArgs(trailingOnly=TRUE)
seed <- as.integer(args[1]) # the seed used for simulation
nsample <- as.integer(args[2]) # the number of sample per condition
method <- as.character(args[3]) # the name of evaluation method
nm <- as.character(args[4]) # the name of used dataset
work_dir <- as.character(args[5]) # the dir containing all files should not end with "/"
work_dir <- paste0(work_dir, "/")

# set R library
.libPaths(paste0(work_dir, "3.5"))

suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(scDD))
suppressPackageStartupMessages(library(MultiAssayExperiment))
suppressPackageStartupMessages(library(survey))

# name of simulated data
part1 <- paste0(nm, "_", nsample, "_", seed, "_part1.rds")
part2 <- paste0(nm, "_", nsample, "_", seed, "_part2.rds")
truth <- paste0(nm, "_", nsample, "_", seed, "_truth.rds")
dataname <- paste0(nm, "_", nsample, "_", seed, ".rds")

cat(paste0(method, ' on ', dataname, '\n'))

datapath <- paste0(work_dir, "simu_data/") # location of simulated data
scriptpath <- paste0(work_dir, "MasterProject/DE_methods/") # location of all DE scripts
resultpath <- paste0(work_dir, "scDD_simu_result/") # location of DE result

# load simulated data
data1 <- readRDS(paste0(datapath, part1))
data2 <- readRDS(paste0(datapath, part2))
condt <- readRDS(paste0(datapath, truth))
count <- rbind(data1$count, data2$count)
tpm <- rbind(data1$tpm, data2$tpm)
data <- list(count=count, tpm=tpm, condt=condt)

# source DE script
source(paste0(scriptpath, method, ".R"))
time <- system.time({res <- get(paste0("run_", method))(data)})

# re-order genes
res$df <- res$df[rownames(data$count),] # re-order genes in result using the order in count matrix

# obtain p and adjusted p
gene_type <- gsub("\\d", "", rownames(res$df))
res$df$gene <- gene_type # gene type
res$df$size <- nsample # sample size

res$time <- time[3] # used time
res$session_info <- sessionInfo() # session info

# save result
saveRDS(res, file = paste0(resultpath, method, "/", method, "_", dataname))

