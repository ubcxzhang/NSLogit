cat("R is running, please do not quit \n")
args = commandArgs(trailingOnly=TRUE)
work_dir <- as.character(args[1]) # the dir containing all files should not end with "/"
work_dir <- paste0(work_dir, "/")

# set R library
.libPaths(paste0(work_dir, "3.5"))

suppressPackageStartupMessages(library(biomaRt))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(AnnotationDbi))
suppressPackageStartupMessages(library(EnsDb.Hsapiens.v75))

edb <- EnsDb.Hsapiens.v75
edb_k <- keys(edb, keytype="SYMBOL")
edb_id <- keys(edb, keytype="GENEID")

# work_dir <- "/home/xsslnc/projects/def-ubcxzh/xsslnc/"
load(paste0(work_dir, "IPFCelltypes.rda")) #IPFList

# load full gene list (include all 0 genes with respect for certain cell type)
G136_gene_dict <- read.delim(paste0(work_dir, "GSE136831_file/GSE136831_AllCells.GeneIDs.txt")) # Ensembl_GeneID, HGNC_EnsemblAlt_GeneID
rownames(G136_gene_dict) <- G136_gene_dict$HGNC_EnsemblAlt_GeneID
G135_gene_SYMBOL <- readRDS(paste0(work_dir, "GSE135893_file/GSE135893_genes.rds"))
G135_gene_SYMBOL <- as.data.frame(G135_gene_SYMBOL)
rownames(G135_gene_SYMBOL) <- G135_gene_SYMBOL[,1]


#G135 <- read.delim(paste0(work_dir, "GSE135893_file/GSE135893_genes.tsv"), header=FALSE)

ensembl <- useEnsembl(
  biomart = "ENSEMBL_MART_ENSEMBL",
  dataset = 'hsapiens_gene_ensembl',
  version = 84)

manual <- as.data.frame(matrix(c(
  c("P11-14N7.2", "ENSG00000232527"),
  c("RP11-123O22.1", "ENSG00000250064"),
  c("CTC-250I14.6", "ENSG00000267598"),
  c("SELM", "ENSG00000198832"),
  c("AC058791.1", "ENSG00000273319"),
  c("RP11-386I14.4", "ENSG00000273338")
), ncol = 2, byrow=TRUE))

colnames(manual) <- c("hgnc_symbol", "ensembl_gene_id")

symb_to_ensg <- function(symb, manual) {
  # use biomart and annotationdbi to map gene symbol to id
  ens <- getBM(
    filters = 'hgnc_symbol',
    values = symb,
    mart = ensembl,
    attributes = c('ensembl_gene_id','hgnc_symbol'))
  mapped <- symb %in% ens$hgnc_symbol
  unfound <- symb[!mapped]
  dbi <- as.character(unfound[unfound %in% edb_k == TRUE])
  if (length(dbi) != 0) {
    dbi_query <- as.data.frame(
      select(edb, keys = dbi, columns=c("GENEID", "SYMBOL"), keytype = "SYMBOL"))
    colnames(dbi_query) <- c("ensembl_gene_id", "hgnc_symbol")
    ens <- rbind(ens, dbi_query)
    mapped <- symb %in% ens$hgnc_symbol
  }

  # add manual mapping
  for (i in 1:dim(manual)[1]) {
    if (manual$hgnc_symbol[i] %in% symb[!mapped]) {
      ens <- rbind(ens, manual[i, c("ensembl_gene_id", "hgnc_symbol")])
    }
  }

  mapped <- symb %in% ens$hgnc_symbol

  # add unmapped gene
  gene_unfound <- as.data.frame(
    cbind(symb[!mapped], rep("", length(symb[!mapped]))))
  colnames(gene_unfound) <- c("hgnc_symbol", "ensembl_gene_id")
  gene <- rbind(ens, gene_unfound)
  mapped <- symb %in% gene$hgnc_symbol
  cat(paste0("contains all gene symbols including those cannot get an id: ", sum(!mapped) == 0), "\n")
  return(gene)
}

ensg_to_symb <- function(ensg, manual) {
  ens <- getBM(
    filters = 'ensembl_gene_id',
    values = ensg,
    mart = ensembl,
    attributes = c('ensembl_gene_id','hgnc_symbol'))
  mapped <- ensg %in% ens$ensembl_gene_id
  unfound <- ensg[!mapped]
  dbi <- as.character(unfound[unfound %in% edb_id == TRUE])
  if (length(dbi) != 0) {
    dbi_query <- as.data.frame(
      select(edb, keys = dbi, columns=c("GENEID", "SYMBOL"), keytype = "GENEID"))
    colnames(dbi_query) <- c("ensembl_gene_id", "hgnc_symbol")
    ens <- rbind(ens, dbi_query)
    mapped <- ensg %in% ens$ensembl_gene_id
  }

  # add manual mapping
  for (i in 1:dim(manual)[1]) {
    if (manual$ensembl_gene_id[i] %in% ensg[!mapped]) {
      ens <- rbind(ens, manual[i, c("ensembl_gene_id", "hgnc_symbol")])
    }
  }

  mapped <- ensg %in% ens$ensembl_gene_id

  # some mapped ensg has na value, set it to ""
  ens[is.na(ens$hgnc_symbol), "hgnc_symbol"] <- ""

  # add unmapped gene
  gene_unfound <- as.data.frame(
    cbind(ensg[!mapped], rep("", length(ensg[!mapped]))))
  colnames(gene_unfound) <- c("ensembl_gene_id", "hgnc_symbol")
  gene <- rbind(ens, gene_unfound)
  mapped <- ensg %in% gene$ensembl_gene_id
  cat(paste0("contains all ensembl id including those cannot get an symbols: ", sum(!mapped) == 0), "\n")
  return(gene)
}

create_136_dict <- function(G136_gene_dict) {
  is_ENSG <- grepl("ENSG", G136_gene_dict$HGNC_EnsemblAlt_GeneID, fixed = TRUE)
  ENSG <- as.character(G136_gene_dict$HGNC_EnsemblAlt_GeneID[is_ENSG])
  SYMB <- as.character(G136_gene_dict$HGNC_EnsemblAlt_GeneID[!is_ENSG])
  map_ENSG <- ensg_to_symb(ENSG, manual) %>% distinct()
  map_SYMB <- G136_gene_dict[SYMB, ]
  colnames(map_SYMB) <- c("ensembl_gene_id", "hgnc_symbol")
  GENE <- rbind(map_SYMB, map_ENSG) %>% distinct()
  GENE$ensembl_gene_id <- as.character(GENE$ensembl_gene_id)
  GENE$hgnc_symbol <- as.character(GENE$hgnc_symbol)
  return(GENE)
}

create_135_dict <- function(G135_gene_SYMBOL) {
  is_ENSG <- grepl("ENSG", G135_gene_SYMBOL[,1], fixed = TRUE) # 0, all are Gene symbol
  GENE <- symb_to_ensg(as.character(G135_gene_SYMBOL[,1]), manual) %>% distinct()
  GENE$ensembl_gene_id <- as.character(GENE$ensembl_gene_id)
  GENE$hgnc_symbol <- as.character(GENE$hgnc_symbol)
  return(GENE)
}

combine_ns_paper_res <- function(data, cell, manual, gene_full_dict, gene_raw_list) {
  # load ns result to get full list of genes under this type
  res <- readRDS(paste0(work_dir, data, "_file/NsLogit_count_result/", cell, ".rds"))

  GENE <- gene_full_dict

  check <- rownames(res$df) %in% GENE$ensembl_gene_id|rownames(res$df) %in% GENE$hgnc_symbol
  cat(paste0("Dict include all ns res gene: ", all(check==TRUE), "\n"))

  GENE$ns_padj <- rep(Inf, dim(GENE)[1])

  # naive loop to fill padj since there may have duplicated gene symbol row (mapped to different gene id) and vice versa
  # a symbol may have multiple ensemble id or a ensemble id may have multiple symbol
  for (i in 1:dim(res$df)[1]) {
    if (rownames(res$df)[i] %in% GENE$ensembl_gene_id == TRUE) {
      GENE[GENE$ensembl_gene_id == rownames(res$df)[i], "ns_padj"] <- res$df[i, "padj"]
    } else {
      GENE[GENE$hgnc_symbol == rownames(res$df)[i], "ns_padj"] <- res$df[i, "padj"]
    }
  }

  # if gene has been filtered out by ns (==0), the ns_padj filed is Inf

  GENE[is.na(GENE$ns_padj), "ns_padj"] <- 1

  check <- rownames(res$df) %in% GENE$ensembl_gene_id|rownames(res$df) %in% GENE$hgnc_symbol
  cat(paste0("Include all ns res gene: ", all(check==TRUE), "\n"))

  if (data == "GSE136831") {
    GENE$paper_wilcox.p <- rep(Inf, dim(GENE)[1])
    GENE$paper_wilcox.fdr <- rep(Inf, dim(GENE)[1])
    GENE$paper_wilcox.bonferroni <- rep(Inf, dim(GENE)[1])

    paper_res <- read.delim(paste0(work_dir, data, "_file/paper_result/media-8.txt"))
    paper_res$gene <- as.character(paper_res$gene)
    check <- paper_res$gene %in% GENE$ensembl_gene_id|paper_res$gene %in% GENE$hgnc_symbol
    cat(paste0("Include all paper res gene: ", all(check==TRUE), "\n"))
    if (all(check==TRUE) == FALSE) {print(paper_res[!check, ])}

    for (i in 1:dim(paper_res)[1]) {
      if (paper_res$gene[i] %in% GENE$ensembl_gene_id == TRUE) {
        GENE[GENE$ensembl_gene_id == paper_res$gene[i], c("paper_wilcox.p", "paper_wilcox.fdr", "paper_wilcox.bonferroni")] <- paper_res[i, c("wilcox.p", "wilcox.fdr", "wilcox.bonferroni")]
      } else {
        GENE[GENE$hgnc_symbol == paper_res$gene[i], c("paper_wilcox.p", "paper_wilcox.fdr", "paper_wilcox.bonferroni")] <- paper_res[i, c("wilcox.p", "wilcox.fdr", "wilcox.bonferroni")]
      }
    }

    GENE$ns.de <- as.integer(GENE$ns_padj < 0.05)
    GENE$paper_wilcox.p.de <- as.integer(GENE$paper_wilcox.p < 0.05)
    GENE$paper_wilcox.fdr.de <- as.integer(GENE$paper_wilcox.fdr < 0.05)
    GENE$paper_wilcox.bonferroni.de <- as.integer(GENE$paper_wilcox.bonferroni < 0.05)

    # indicate filtered gene
    GENE[GENE$ns_padj == Inf, "ns.de"] <- -1
    GENE[GENE$paper_wilcox.p == Inf, "paper_wilcox.p.de"] <- -1
    GENE[GENE$paper_wilcox.fdr == Inf, "paper_wilcox.fdr.de"] <- -1
    GENE[GENE$paper_wilcox.bonferroni == Inf, "paper_wilcox.bonferroni.de"] <- -1

  } else {

    GENE$paper_p_val <- rep(Inf, dim(GENE)[1])
    GENE$paper_p_val_adj <- rep(Inf, dim(GENE)[1])
    paper_res <- read.csv(file = paste0(work_dir, "GSE135893_file/Disease_vs_Control/", cell, "_disease_vs_control_.csv"),  header = TRUE)

    check <- paper_res$X %in% GENE$ensembl_gene_id|paper_res$X %in% GENE$hgnc_symbol
    cat(paste0("Include all paper res gene: ", all(check==TRUE), "\n"))
    if (all(check==TRUE) == FALSE) {print(paper_res[!check, ])}

    for (i in 1:dim(paper_res)[1]) {
      if (paper_res$X[i] %in% GENE$ensembl_gene_id == TRUE) {
        GENE[GENE$ensembl_gene_id == paper_res$X[i], c("paper_p_val", "paper_p_val_adj")] <- paper_res[i, c("p_val", "p_val_adj")]
      } else {
        GENE[GENE$hgnc_symbol == paper_res$X[i], c("paper_p_val", "paper_p_val_adj")] <- paper_res[i, c("p_val", "p_val_adj")]
      }
    }

    GENE$ns.de <- as.integer(GENE$ns_padj < 0.05)
    GENE$paper_p_val.de <- as.integer(GENE$paper_p_val < 0.05)
    GENE$paper_p_val_adj.de <- as.integer(GENE$paper_p_val_adj < 0.05)

    # indicate filtered gene
    GENE[GENE$ns_padj == Inf, "ns.de"] <- -1
    GENE[GENE$paper_p_val == Inf, "paper_p_val.de"] <- -1
    GENE[GENE$paper_p_val_adj == Inf, "paper_p_val_adj.de"] <- -1
  }

  GENE$cell_type <- cell
  GENE$dataset <- data

  return(GENE)
}


# main

# create mapping using full gene list
gene_136_dict <- create_136_dict(G136_gene_dict)
gene_135_dict <- create_135_dict(G135_gene_SYMBOL)

for (i in 1:length(IPFList$GSE136831)) {
  print(IPFList$GSE136831[i])
  GSE136831_cell <- IPFList$GSE136831[i]
  GSE136831_res <- combine_ns_paper_res("GSE136831", GSE136831_cell, manual, gene_136_dict, G136_gene_dict)
  GSE135893_cell <- gsub("\\/| ", "_", IPFList$GSE135893[i])
  GSE135893_res <- combine_ns_paper_res("GSE135893", GSE135893_cell, manual, gene_135_dict, G135_gene_SYMBOL)

  # add column to indicate whether this rwo appears in another files by id/symbol
  GSE136831_res$id_in_135 <- GSE136831_res$ensembl_gene_id %in% GSE135893_res$ensembl_gene_id
  GSE136831_res[GSE136831_res$ensembl_gene_id == "", "id_in_135"] <- FALSE
  GSE136831_res$sb_in_135 <- GSE136831_res$hgnc_symbol %in% GSE135893_res$hgnc_symbol
  GSE136831_res[GSE136831_res$hgnc_symbol == "", "sb_in_135"] <- FALSE
  GSE136831_res$gene_in_135 <- GSE136831_res$id_in_135|GSE136831_res$sb_in_135

  GSE135893_res$id_in_136 <- GSE135893_res$ensembl_gene_id %in% GSE136831_res$ensembl_gene_id
  GSE135893_res[GSE135893_res$ensembl_gene_id == "", "id_in_136"] <- FALSE
  GSE135893_res$sb_in_136 <- GSE135893_res$hgnc_symbol %in% GSE136831_res$hgnc_symbol
  GSE135893_res[GSE135893_res$hgnc_symbol == "", "sb_in_136"] <- FALSE
  GSE135893_res$gene_in_136 <- GSE135893_res$id_in_136|GSE135893_res$sb_in_136

  # add column for de in the other dataset; -2 if this gene not exisits
  GSE135893_res$ns136.de <- rep(-2, dim(GSE135893_res)[1])
  GSE135893_res$paper136_wilcox.p.de <- rep(-2, dim(GSE135893_res)[1])
  GSE135893_res$paper136_wilcox.fdr.de <- rep(-2, dim(GSE135893_res)[1])
  GSE135893_res$paper136_wilcox.bonferroni.de <- rep(-2, dim(GSE135893_res)[1])

  for (i in 1:dim(GSE135893_res)[1]) {
    # -2 for unfound gene
    if (GSE135893_res[i, "gene_in_136"] == FALSE) { next() }
    if (GSE135893_res[i, "id_in_136"] == TRUE) { # by ensemble id
      find_136 <- GSE136831_res %>% dplyr::filter(ensembl_gene_id == GSE135893_res[i, "ensembl_gene_id"])
    } else { # by gene symbol
      find_136 <- GSE136831_res %>% dplyr::filter(hgnc_symbol == GSE135893_res[i, "hgnc_symbol"])
    }
    GSE135893_res[i, c("ns136.de", "paper136_wilcox.p.de", "paper136_wilcox.fdr.de", "paper136_wilcox.bonferroni.de")] <-
      find_136[1, c("ns.de", "paper_wilcox.p.de", "paper_wilcox.fdr.de", "paper_wilcox.bonferroni.de")]
  }

  GSE136831_res$ns135.de <- rep(-2, dim(GSE136831_res)[1])
  GSE136831_res$paper135_p_val.de <- rep(-2, dim(GSE136831_res)[1])
  GSE136831_res$paper135_p_val_adj.de <- rep(-2, dim(GSE136831_res)[1])

  for (i in 1:dim(GSE136831_res)[1]) {
    # -2 for unfound gene
    if (GSE136831_res[i, "gene_in_135"] == FALSE) { next() }
    if (GSE136831_res[i, "id_in_135"] == TRUE) { # by ensemble id
      find_135 <- GSE135893_res %>% dplyr::filter(ensembl_gene_id == GSE136831_res[i, "ensembl_gene_id"])
    } else { # by gene symbol
      find_135 <- GSE135893_res %>% dplyr::filter(hgnc_symbol == GSE136831_res[i, "hgnc_symbol"])
    }
    GSE136831_res[i, c("ns135.de", "paper135_p_val.de", "paper135_p_val_adj.de")] <-
      find_135[1, c("ns.de", "paper_p_val.de", "paper_p_val_adj.de")]
  }

  write.csv(GSE136831_res, file = paste0(work_dir, "real_data_comb_res/GSE136831/GSE136831_", GSE136831_cell, ".csv"))
  write.csv(GSE135893_res, file = paste0(work_dir, "real_data_comb_res/GSE135893/GSE135893_", GSE135893_cell, ".csv"))
}
