cat("R is running, please do not quit \n")
args = commandArgs(trailingOnly=TRUE)
work_dir <- as.character(args[1]) # the dir containing all files should not end with "/"
work_dir <- paste0(work_dir, "/")

# set R library
# .libPaths(paste0(work_dir, "3.5"))

suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))

realdata_file <- paste0(work_dir, "GSE135893_file/")
meta <- read.csv(paste0(realdata_file, "GSE135893_IPF_metadata.csv"), header = TRUE)
ns_counts_path <- paste0(work_dir, "GSE135893_file/NsLogit_count_result/")
profile <- paste0(work_dir, "GSE135893_count/")
gene_list <- paste0(work_dir, "gene_list/")

paper_file <- paste0(work_dir, "GSE135893_file/Disease_vs_Control/")
cells <- as.character(unique(meta$celltype))
compare <- data.frame(); hist_data <- data.frame()

all_files <- list.files(path = ns_counts_path, pattern = "rds$")

i <- 1
for (c in cells) {
  if (c == "HAS1 High Fibroblasts") {next()}
  files <- gsub("\\/| ", "_", c)
  ns_file <- paste0(files, ".rds")
  if (ns_file %in% all_files == FALSE) {next()}

  print(files)
  ns_res <- readRDS(paste0(ns_counts_path, files, ".rds")) # ns result
  se_res <- read.csv(file = paste0(paper_file, files, "_disease_vs_control_.csv"),  header = TRUE) # paper result

  rownames(se_res) <- se_res$X
  ns_res <- ns_res$df
  ns_res$gene <- rownames(ns_res)

  gene_profile <- readRDS(paste0(profile, files, "_fckl.rds"))
  ns_res <- cbind(ns_res, gene_profile)

  # calculate number of gene with logfc > 0.25
  top_range <- sum(ns_res$fc > 0.25)

  # order result with descending kl_diff
  ns_res <- ns_res[order(-ns_res$kl_diff),]

  # top significant gene for kl_diff
  top_gene <- rownames(ns_res[1:top_range,])

  # filter NS result
  valid_1 <- ns_res %>% filter(fc > 0.25)
  valid_2 <- ns_res %>% filter(fc < 0.25&gene %in% top_gene)
  select_gene <- c(valid_1$gene, valid_2$gene)
  ns_res <- ns_res[select_gene, ]

  # re-calculate
  ns_res$padj <- p.adjust(ns_res$pval, method="fdr", n=nrow(ns_res))

  # compare
  ns_res[is.na(ns_res$padj), "padj"] <- 1
  se_res[is.na(se_res$p_val_adj), "p_val_adj"] <- 1
  ns_gene <- rownames(ns_res[ns_res$padj<0.05,])
  se_gene <- rownames(se_res[se_res$p_val_adj<0.05,])

  miss <- vector()
  novel <- vector()
  novel_rep <- vector()
  for (x in ns_gene) {if (x %in% se_gene == FALSE) {novel <- c(novel, x)}}
  for (x in se_gene) {if (x %in% ns_gene == FALSE) {miss <- c(miss, x)}}

  # novel and not reported by paper: due to filer
  for (x in novel) {if (x %in% rownames(se_res) == FALSE) {novel_rep <- c(novel_rep, x)}}

  compare[c, "NsLogit"] <- length(ns_gene)
  compare[c, "Published"] <- length(se_gene)
  compare[c, "NsLogit_Novel"] <- length(novel)
  compare[c, "NsLogit_Miss"] <- length(miss)
  compare[c, "Novel_from_filter"] <- length(novel_rep)
  compare[c, "take top(%) of kl_diff"] <- top_range/dim(gene_profile)[1]
  i <- i+2

  colnames(ns_res) <- paste0("NsLogit_", colnames(ns_res))
  colnames(se_res) <- paste0("Published_", colnames(se_res))

  ns_list <- ns_res[ns_gene, ]
  se_list <- se_res[se_gene, ]

  write.csv(se_list, file = paste0(gene_list, files, "/", files, "_Published.csv"))
  write.csv(ns_list, file = paste0(gene_list, files, "/", files, "_NsLogit.csv"))

  empty_res <- matrix(NA, nrow = length(novel_rep), ncol = dim(se_res)[2])
  rownames(empty_res) <- novel_rep
  colnames(empty_res) <- colnames(se_res)
  empty_res <- as.data.frame(empty_res)
  se_res <- rbind(se_res, empty_res)

  novel_list <- cbind(ns_res[novel, ], se_res[novel, ])
  novel_list$not_in_paper <- rownames(novel_list) %in% novel_rep

  write.csv(novel_list, file = paste0(gene_list, files, "/", files, "_NsLogit_Novel.csv"))

  # plot not in paper novel
  data <- readRDS(paste0(profile, files, ".rds"))
  for (g in novel_rep) {
    plot_name <- paste0(g, "_", files)
    g_data <- data.frame(count=log(data$count[g, ] + 1), group=data$condt)
    g_data <- g_data[g_data$count != 0,]
    p <- ggplot(g_data, aes(x=count, color=group)) +
      geom_density() +
      labs(title = paste0("Gene: ", g, " Cell: ", files), x = "non-zero log(count + 1)")
    png(filename = paste0(gene_list, "plot/", plot_name, ".png"), width = 5, height = 5, units = 'in', res = 300)
    grid.arrange(p, ncol = 1)
    dev.off()
  }
}

print(compare)
write.csv(compare, file = paste0(gene_list, "compare.csv"))


