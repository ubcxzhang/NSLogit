args = commandArgs(trailingOnly=TRUE)
filename <- as.character(args[1]) # the cell subtype
work_dir <- as.character(args[2]) # the dir containing all files should not end with "/"
work_dir <- paste0(work_dir, "/")

# set R library
.libPaths(paste0(work_dir, "3.5"))
library(Matrix)

# original file of GSE136831
# downloadable from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE136831

realdata_file <- paste0(work_dir, "GSE136831_file/")
realdata_part_file <- paste0(work_dir, "GSE136831_count/")

# load full data
data <- readMM(paste0(realdata_file, "GSE136831_RawCounts_Sparse.mtx")) # [1]  45947 312928
cell <- read.delim(paste0(realdata_file,"GSE136831_AllCells.cellBarcodes.txt"), header = FALSE) # [1] 312928      1
gene <- read.delim(paste0(realdata_file,"GSE136831_AllCells.GeneIDs.txt")) # [1] 45947     2
meta <- read.delim(paste0(realdata_file,"GSE136831_AllCells.Samples.CellType.MetadataTable.txt")) # [1] 312928      9
#data <- as.data.frame(as.matrix(data))
rownames(data) <- gene$HGNC_EnsemblAlt_GeneID
colnames(data) <- cell[,1]
head(data)
# select cell
select <- meta[meta$Manuscript_Identity == filename, ]
dim(select)
head(select)
condt <- select$Disease_Identity
head(condt)
select_cell <- select$CellBarcode_Identity
head(select_cell)
count <- data[ ,select_cell]

# remove all-zero genes
count <- count[rowSums(count) > 0, ]
sct <- list(count=count, condt=condt)

saveRDS(sct, file = paste0(realdata_part_file, filename, ".rds"))
