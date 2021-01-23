cat("R is running, please do not quit \n")
args = commandArgs(trailingOnly=TRUE)
work_dir <- as.character(args[1]) # the dir containing all files should not end with "/"
work_dir <- paste0(work_dir, "/")

# set R library
.libPaths(paste0(work_dir, "3.5"))

suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))

realdata_file <- paste0(work_dir, "GSE135893_file/")
meta <- read.csv(paste0(realdata_file, "GSE135893_IPF_metadata.csv"), header = TRUE)
ns_counts_path <- paste0(work_dir, "GSE135893_file/NsLogit_count_result/")

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
  
  ns_res <- readRDS(paste0(ns_counts_path, files, ".rds")) # ns result
  se_res <- read.csv(file = paste0(paper_file, files, "_disease_vs_control_.csv"),  header = TRUE) # paper result

  rownames(se_res) <- se_res$X

  # paper result only report genes that pass their filter criteria: see Seurat_NB.R
  # subset ns result == apply the same filter

  ns_res <- ns_res$df[rownames(se_res), ]
  ns_res[is.na(ns_res$padj), "padj"] <- 1
  se_res[is.na(se_res$p_val_adj), "p_val_adj"] <- 1
  ns_gene <- rownames(ns_res[ns_res$padj<0.05,])
  se_gene <- rownames(se_res[se_res$p_val_adj<0.05,])

  miss <- vector()
  novel <- vector()
  for (x in ns_gene) {if (x %in% se_gene == FALSE) {novel <- c(novel, x)}}
  for (x in se_gene) {if (x %in% ns_gene == FALSE) {miss <- c(miss, x)}}

  compare[c, "NsLogit"] <- length(ns_gene)
  compare[c, "Published"] <- length(se_gene)
  compare[c, "NsLogit_Novel"] <- length(novel)
  compare[c, "NsLogit_Miss"] <- length(miss)
  hist_data[i, "Cell"] <- c
  hist_data[i, "DE"] <- length(ns_gene)
  hist_data[i, "Method"] <- "NsLogit"
  hist_data[(i+1), "Cell"] <- c
  hist_data[(i+1), "DE"] <- length(se_gene)
  hist_data[(i+1), "Method"] <- "Published"
  i <- i+2
}

print(compare)

pop <- cbind(as.character(meta$celltype), as.character(meta$population))
pop <- unique(pop); rownames(pop) <- pop[,1]
hist_data$population <- pop[hist_data$Cell, 2]

plot_stat <- function(data, text) {
  dodge <- position_dodge(width=0.9)
  p <- ggplot(data, aes(x = reorder(Cell,-DE), y = DE), fill=as.factor(Method)) +
    geom_bar(stat="identity", color="white", position=dodge, aes(fill = Method)) +
    scale_fill_manual(name="Method", values=c('#373F47','#8B8982')) +
    labs(title=text, x="Cell Types", y = "Number of DE Genes")+
    facet_grid(~population,
      scales = "free",
      space = "free")+
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
      legend.position="bottom", legend.box = "horizontal",
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14),
      strip.text.x = element_text(size = 9),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black"))
  p
}

p <- plot_stat(hist_data, "")
png(filename = paste0(work_dir, "Figure/real.png"), width = 10, height = 5, units = 'in', res = 100)
grid.arrange(p, ncol = 1)
dev.off()

cat("Figure real.png saved\n")
cat("Done!\n")

