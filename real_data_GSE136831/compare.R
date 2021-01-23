cat("R is running, please do not quit \n")
args = commandArgs(trailingOnly=TRUE)
work_dir <- as.character(args[1]) # the dir containing all files should not end with "/"
work_dir <- paste0(work_dir, "/")

# work_dir <- "/home/xsslnc/projects/def-ubcxzh/xsslnc/"
load(paste0(work_dir, "IPFCelltypes.rda")) #IPFList
library(Matrix)
res <- list()

.libPaths(paste0(work_dir, "3.5"))

suppressPackageStartupMessages(library(biomaRt))
suppressPackageStartupMessages(library(AnnotationDbi))
suppressPackageStartupMessages(library(EnsDb.Hsapiens.v75))

edb <- EnsDb.Hsapiens.v75
edb_k <- keys(edb, keytype="SYMBOL")

manual <- as.data.frame(matrix(c(
  c("P11-14N7.2", "ENSG00000232527"),
  c("RP11-123O22.1", "ENSG00000250064"),
  c("CTC-250I14.6", "ENSG00000267598"),
  c("SELM", "ENSG00000198832"),
  c("AC058791.1", "ENSG00000273319"),
  c("RP11-386I14.4", "ENSG00000273338")
), ncol = 2, byrow=TRUE))

colnames(manual) <- c("hgnc_symbol", "ensembl_gene_id")

realdata_part_file <- paste0(work_dir, "GSE135893_count/")
sct <- readRDS(paste0(realdata_part_file, filename, ".rds"))

id <- rownames(sct$count)

ensembl <- useEnsembl(
  biomart = "ENSEMBL_MART_ENSEMBL",
  dataset = 'hsapiens_gene_ensembl',
  version = 84)


compare <- data.frame()

for (i in 1:length(IPFList$GSE136831)) {
  cur <- data.frame()

  print(IPFList$GSE136831[i])
  GSE136831_cell <- IPFList$GSE136831[i]
  GSE135893_cell <- files <- gsub("\\/| ", "_", IPFList$GSE135893[i])
  GSE135893_novel <- read.csv(
    file = paste0(
      work_dir, "gene_list/", GSE135893_cell, "/",
      GSE135893_cell, "_NsLogit_Novel.csv"),
    header = TRUE)

  if (nrow(GSE135893_novel) == 0) {next()}
  GSE136831 <- readRDS(paste0(work_dir, "GSE136831_count/", GSE136831_cell, ".rds"))

  # GSE136831 use both ensemble id and gene symbol
  is_ENSG <- grepl("ENSG", rownames(GSE136831$count), fixed = TRUE)
  Novel_Gene <- GSE135893_novel$NsLogit_gene
  GSE136831_ENSG <- rownames(GSE136831$count)[is_ENSG]
  GSE136831_SYMB <- rownames(GSE136831$count)[!is_ENSG]

  GSE135893_novel_ens <- getBM(
    filters = 'hgnc_symbol',
    values = Novel_Gene,
    mart = ensembl,
    attributes = c('ensembl_gene_id','hgnc_symbol'))
  GSE135893_novel_ens <- rbind(GSE135893_novel_ens, manual)

  novel_has_ens <- Novel_Gene %in% GSE135893_novel_ens$hgnc_symbol

  if (!(all(novel_has_ens) == TRUE)) {
    # use AnnotationDbi to find unfound ensemble id of gene symbol
    novel_unfound <- Novel_Gene[!novel_has_ens]
    k <- as.character(novel_unfound[novel_unfound %in% edb_k == TRUE])
    if (length(k) != 0) {
      query <- as.data.frame(
        select(edb, keys = k, columns=c("GENEID", "SYMBOL"), keytype = "SYMBOL"))
      colnames(query) <- c("ensembl_gene_id", "hgnc_symbol")
      GSE135893_novel_ens <- rbind(GSE135893_novel_ens, query)

      novel_has_ens <- Novel_Gene %in% GSE135893_novel_ens$hgnc_symbol
    }
  }

  if (!(all(novel_has_ens) == TRUE)) {
    cat("following novel genes has not been mapped to ensembl id: ")
    print(as.character(Novel_Gene[!novel_has_ens]))
  } else {
    cat("All novel genes has been mapped to ensembl id \n")
  }

  GSE136831_biomart <- getBM(
    filters = 'hgnc_symbol',
    values = GSE136831_SYMB,
    mart = ensembl,
    attributes = c('ensembl_gene_id','hgnc_symbol'))

  # in case that not all symbol from GSE136831 find their emsemble id
  gse136_has_ens <- GSE136831_SYMB %in% GSE136831_biomart$hgnc_symbol

  if (!(all(gse136_has_ens) == TRUE)) {
    # use AnnotationDbi to find unfound ensemble id of gene symbol
    gse136_unfound <- GSE136831_SYMB[!gse136_has_ens]
    k <- as.character(gse136_unfound[gse136_unfound %in% edb_k == TRUE])
    if (length(k) != 0) {
      query <- as.data.frame(
        select(edb, keys = k, columns=c("GENEID", "SYMBOL"), keytype = "SYMBOL"))
      colnames(query) <- c("ensembl_gene_id", "hgnc_symbol")
      GSE136831_biomart <- rbind(GSE136831_biomart, query)
      gse136_has_ens <- GSE136831_SYMB %in% GSE136831_biomart$hgnc_symbol
    }
  }

  # keep gene for current type (since manual include other types)
  GSE135893_novel_ens <- GSE135893_novel_ens[GSE135893_novel_ens$hgnc_symbol %in% Novel_Gene,]

  ens <- GSE135893_novel_ens$ensembl_gene_id %in% GSE136831_ENSG
  mapped_ens <- GSE135893_novel_ens$ensembl_gene_id %in% GSE136831_biomart$ensembl_gene_id

  valid <- 0

  if (all(ens) == FALSE) {
    cat("No overlap between ensembl id of GSE135893_novel and GSE136831's ensembl id \n")
  } else {
    print(as.character(Novel_Gene[!is.na(Novel_Gene[ens])]))
    valid <- valid + length(Novel_Gene[!is.na(Novel_Gene[ens])])
  }

  if (all(mapped_ens) == FALSE) {
    cat("No overlap between ensembl id of GSE135893_novel and GSE136831's mapped ensembl id \n")
  } else {
    print(as.character(Novel_Gene[!is.na(Novel_Gene[mapped_ens])]))
    valid <- valid + length(Novel_Gene[!is.na(Novel_Gene[mapped_ens])])
  }

  unmapped_gse136 <- GSE136831_SYMB[!gse136_has_ens]

  cur[1, "Type"] <- IPFList$GSE136831[i]
  cur[1, "#Novel Gene"] <- length(Novel_Gene)
  cur[1, "#Valided by GSE136831"] <- valid
  cur[1, "#unmapped GSE136831 Gene"] <- length(unmapped_gse136)
  compare <- rbind(compare, cur)
}
