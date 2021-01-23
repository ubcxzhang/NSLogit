args = commandArgs(trailingOnly=TRUE)
filename <- as.character(args[1])
work_dir <- as.character(args[2]) # the dir containing all files should not end with "/"
work_dir <- paste0(work_dir, "/")

# set R library
.libPaths(paste0(work_dir, "3.5"))

# process GSE135893_ILD_annotated_fullsize.rds into subfiles
# each contains a list for a cell type
# need Seurat 3.0.0 to run
library(Seurat)
library(dplyr)

# original file of GSE135893
# downloadable from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE135893

realdata_file <- paste0(work_dir, "GSE135893_file/")
realdata_part_file <- paste0(work_dir, "GSE135893_count/")

meta <- read.csv(paste0(realdata_file, "GSE135893_IPF_metadata.csv"), header = TRUE)
cellname <- as.character(unique(meta$celltype))
names(cellname) <- gsub(" |/", "_", cellname)

# load fullsize file
ild <- readRDS(paste0(realdata_file, "GSE135893_ILD_annotated_fullsize.rds"))
i <- cellname[filename] # get real cell name from its legal file name
# subset a certain cell type
cells <- SubsetData(ild, cells = row.names(ild@meta.data[ild@meta.data$celltype == i, ]))
# get its count matrix
data.use <- GetAssayData(object = cells[["SCT"]], slot = "counts")
count <- as.matrix(data.use)
condt <- cells@meta.data$Status
count <- count[rowSums(count) > 0, ]
sct <- list(count=count, condt=condt)

saveRDS(sct, file = paste0(realdata_part_file, filename, ".rds"))

