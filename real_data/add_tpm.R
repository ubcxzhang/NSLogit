args = commandArgs(trailingOnly=TRUE)
filename <- as.character(args[1])
work_dir <- as.character(args[2]) # the dir containing all files should not end with "/"
work_dir <- paste0(work_dir, "/")

# set R library
.libPaths(paste0(work_dir, "3.5"))

suppressPackageStartupMessages(library(biomaRt))
suppressPackageStartupMessages(library(AnnotationDbi))
suppressPackageStartupMessages(library(EnsDb.Hsapiens.v75))

edb <- EnsDb.Hsapiens.v75

realdata_part_file <- paste0(work_dir, "GSE135893_count/")
sct <- readRDS(paste0(realdata_part_file, filename, ".rds"))

id <- rownames(sct$count)

ensembl <- useEnsembl(
  biomart = "ENSEMBL_MART_ENSEMBL",
  dataset = 'hsapiens_gene_ensembl',
  version = 84)
biomart_found <- getBM(
  filters = "hgnc_symbol",
  values = id,
  mart = ensembl,
  attributes = c('start_position', 'end_position','hgnc_symbol'))
biomart_found <- as.data.frame(biomart_found)
biomart_found$length <- biomart_found$end_position - biomart_found$start_position
bm <- aggregate(biomart_found[, c("length")], list(biomart_found$"hgnc_symbol"), max)
colnames(bm) <- c("SYMBOL", "length")
rownames(bm) <- bm$"SYMBOL"

unfound <- id[id %in% bm$SYMBOL == FALSE]

edb_k <- keys(edb, keytype="SYMBOL")
k <- unfound[unfound %in% edb_k == TRUE]
query <- select(edb, keys = k, columns=c("GENESEQEND", "GENESEQSTART", "SYMBOL"),
  keytype = "SYMBOL")
query <- as.data.frame(query)
query$length <- query$GENESEQEND - query$GENESEQSTART
annodbi <- aggregate(query[, c("length")], list(query$"SYMBOL"), max)
colnames(annodbi) <- c("SYMBOL", "length")
rownames(annodbi) <- annodbi$SYMBOL

gene <- rbind(bm, annodbi)

cat(paste0("Cell type ", filename, " has ", length(id), " non-zero genes, ", dim(gene)[1],
  " of them have gene length and are kept for DE. \n"))

paper_file <- paste0(work_dir, "GSE135893_file/Disease_vs_Control/")
se_res <- read.csv(file = paste0(paper_file, filename, "_disease_vs_control_.csv"),  header = TRUE) # paper result

miss <- se_res[se_res$X %in% gene$SYMBOL == FALSE, "X"]

cat(paste0(length(miss), " genes that do not have gene length appear in paper's result. \n"))

sct$count <- sct$count[gene$SYMBOL, ]
cat("Only genes with a gene length are kept in count. \n")
sct$tpm <- sct$count

# normalize for gene length
for (g in gene$SYMBOL) {sct$tpm[g, ] <- sct$count[g, ]/gene[g, "length"]}
# normalize for gene depth
colsum <- colSums(sct$tpm)
for (i in 1:length(colsum)) {sct$tpm[, i] <- sct$tpm[, i]/colsum[i]*1000000}

realdata_part_file <- paste0(work_dir, "GSE135893_with_tpm/")
saveRDS(sct, file = paste0(realdata_part_file, filename, ".rds"))

cat(paste0("TPM is calculated\n",  realdata_part_file, filename, ".rds is saved\n"))

# look into gene with many length
gene_list <- matrix(0, nrow = length(gene$SYMBOL), ncol = 6)
colnames(gene_list) <- c("count", "min", "max", "max-min", "median", "(max-min)/median")
rownames(gene_list) <- gene$SYMBOL
biomart_found$SYMBOL <- biomart_found$hgnc_symbol
biomart_found <- biomart_found[, c("SYMBOL", "length")]
query <- query[, c("SYMBOL", "length")]
fullres <- rbind(biomart_found, query)
for (g in gene$SYMBOL) {
  gene_list[g, "count"] <- sum(fullres$SYMBOL == g)
  if (sum(fullres$SYMBOL == g) == 1) {next()}
  sel_res <- fullres[fullres$SYMBOL == g, ]
  gene_list[g, "min"] <- min(sel_res$length)
  gene_list[g, "max"] <- max(sel_res$length)
  gene_list[g, "max-min"] <- max(sel_res$length) - min(sel_res$length)
  gene_list[g, "median"] <- median(sel_res$length)
}

gene_list <- as.data.frame(gene_list)
gene_list <- gene_list[gene_list$"count" != 1 & gene_list$"max-min" != 0,]
gene_list$"(max-min)/median" <- gene_list$"max-min"/gene_list$"median"
cat("Summary on genes with more than one records with different length (we take the max): \n")
gene_list[order(gene_list$"(max-min)/median"),]


