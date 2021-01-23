
set.seed(12345)
# parameter
args = commandArgs(trailingOnly=TRUE)
dataname <- as.character(args[1]) # the name of original dataset to use
nsample <- as.integer(args[2]) # the number of cell per gene
seed <- as.integer(args[3]) # seed for simulation
work_dir <- as.character(args[4]) # the dir containing all files should not end with "/"
work_dir <- paste0(work_dir, "/")

# set R library
.libPaths(paste0(work_dir, "3.5"))

# configuration
if (dataname == "GSE74596") {
  groupid = "source_name_ch1"
  keepgroup = c("Single_cell_RNA-seq_NKT0", "Single_cell_RNA-seq_NKT17")
} else if (dataname == "GSE45719") {
  groupid = "source_name_ch1"
  keepgroup = c("16-cell stage blastomere", "Mid blastocyst cell (92-94h post-fertilization)")
} else if (dataname == "GSE60749-GPL13112") {
  groupid = c("source_name_ch1", "characteristics_ch1.1")
  keepgroup = c("v6.5 mouse embryonic stem cells.culture conditions: 2i+LIF", "v6.5 mouse embryonic stem cells.culture conditions: serum+LIF")
}

# location of simulation_scDD.R
source(paste0(work_dir, "MasterProject/scDD_framework/simulation_scDD.R"))

# location of original dataset
datapath <- paste0(work_dir, "conquer_data/")

# location of result path
resultpath <- work_dir

# simulation
sce <- format_scDD(datapath, dataname, groupid, keepgroup)
SD <- simulation_scDD(sce, nsample, seed)
condt <- colData(SD)$condition
data <- list(count=assay(SD))
data$count <- round(data$count) # round count to integer
# calculate TPM: assume all genes are equal length
col_div <- function(dat, col_sum) {
  for (i in 1:length(col_sum)) {dat[,i] <- dat[,i]/col_sum[i]}
  dat
}
half <- dim(data$count)[1]/2; full <- dim(data$count)[1]
data$tpm <- col_div(data$count, colSums(data$count))*1000000

data1 <- list(count=data$count[1:half,], tpm=data$tpm[1:half,])
data2 <- list(count=data$count[(half+1):full,], tpm=data$tpm[(half+1):full,])

saveRDS(data1, file = paste0(resultpath, "simu_data/", dataname, "_", nsample, "_", seed, "_part1.rds"))
saveRDS(data2, file = paste0(resultpath, "simu_data/", dataname, "_", nsample, "_", seed, "_part2.rds"))
saveRDS(condt, file = paste0(resultpath, "simu_data/", dataname, "_", nsample, "_", seed, "_truth.rds"))

