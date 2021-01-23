args = commandArgs(trailingOnly=TRUE)
work_dir <- as.character(args[1]) # the dir containing all files should not end with "/"
work_dir <- paste0(work_dir, "/")

# set R library
.libPaths(paste0(work_dir, "3.5"))

source(paste0(work_dir, "MasterProject/utils.R"))

# conquer sim123 data path, and generated subfile will also be saved there
datapath <- paste0(work_dir, "conquer_data/")

data <- read_sim_data_from_rds(datapath, "GSE74596sim123")
gen_subfile(data, 42, c(44, 22, 12, 6), c(1, 5, 5, 5), datapath, "GSE74596sim123_subfile")

data <- read_sim_data_from_rds(datapath, "GSE60749-GPL13112sim123")
gen_subfile(data, 42, c(90, 48, 24, 12), c(1, 5, 5, 5), datapath, "GSE60749-GPL13112sim123_subfile")

data <- read_sim_data_from_rds(datapath, "GSE45719sim123")
gen_subfile(data, 42, c(50, 24, 12, 6), c(1, 5, 5, 5), datapath, "GSE45719sim123_subfile")

