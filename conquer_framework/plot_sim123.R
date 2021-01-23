args = commandArgs(trailingOnly=TRUE)
work_dir <- as.character(args[1]) # the dir containing all files should not end with "/"
work_dir <- paste0(work_dir, "/")

# set R library
.libPaths(paste0(work_dir, "3.5"))

suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(cowplot))

#' @param res result dataframe
#' @param text title text
#' @param rec shaded area information
#' @param show.size show the number of points if not NA; show n by set a y value
#'
plot_all <- function(res, text, rec, show.size) {
  # point color
  point_col <- c("#800E13", "#1F78B4", "#666666", "#B15928",
    "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6",
    "#FBBFCA", "#AACC00", "#2EC4B6", "#739E82", "#2C5530", "#34623F", "#800E13", "#49306B")
  names(point_col) <- c("50", "100", "500", "1000", "1500", "2000", "2500", "3000", "3500",
    "6", "12", "22", "24", "44", "48", "50", "90")

  p <- ggplot()+
    scale_y_discrete(dim(rec)[1]) +
    geom_rect(data = rec, aes(ymin = xmin, ymax = xmax, alpha = alpha),
      xmin = -Inf, xmax = Inf, fill = "#E6E6EA")+
    geom_boxplot(data=res, aes(y = Method, x = Value), outlier.shape = "")+
    geom_jitter(data=res, size=1.5, shape=16, position=position_jitter(0.2), aes(y = Method, x = Value, colour=Size))+
    labs(title=paste0(text), x = "", y = "", size=25)+
    scale_color_manual(values = point_col)+
    #facet_grid(~Parametric,
    #  scales = "free",
    #  space = "free")+
    xlim(0,1) +
    theme(
      axis.text.x = element_text(size=15), # y axis text size
      axis.text.y = element_text(size=20), # x axis text size
      plot.title = element_text(size=25,face="bold"),
      legend.text = element_text(size=15), # legend text size
      panel.background = element_rect(fill = "white", colour = "black"),
      legend.position="bottom", legend.box = "horizontal",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()) +
    guides(alpha = FALSE)

  if (is.na(show.size)==FALSE) {
    p <- p + geom_text(data=res, aes(x = Method, y = show.size, label = n), angle = 90, size = 3.5)
  }

  p
}

rename <- function(res, old, new) {
  res$FDR[res$FDR$Method==old, "Method"] <- new
  res$TPR[res$TPR$Method==old, "Method"] <- new
  res$F1[res$F1$Method==old, "Method"] <- new
  return(res)
}

order_method <- function(res, set_order) {
  res$FDR$Method <- factor(res$FDR$Method,
    levels = set_order, ordered = TRUE)
  res$TPR$Method <- factor(res$TPR$Method,
    levels = set_order, ordered = TRUE)
  res$F1$Method <- factor(res$F1$Method,
    levels = set_order, ordered = TRUE)
  return(res)
}

unfilter_order <- c(
  "NsLogit", "Logistic",
  "SCDE", "scDD", "D3E",
  "ttest", "Wilcoxon", "limmatrend",
  "ROTSvoom", "ROTScpm", "SeuratBimodIsExpr2",
  "edgeRLRTcensus", "DESeq2census", "ROTStpm",
  "SAMseq", "MASTcpm", "MASTcpmDetRate",
  "MASTtpm", "voomlimma", "SeuratTobit",
  "metagenomeSeq", "MASTtpmDetRate", "DEsingle",
  "monoclecensus", "DESeq2nofilt", "monoclecount",
  "edgeRLRT", "edgeRLRTdeconv", "NODES",
  "DESeq2", "BPSC", "DESeq2betapFALSE",
  "edgeRLRTrobust", "edgeRQLFDetRate", "edgeRQLF",
  "SeuratBimod", "monocle", "SeuratBimodnofilt")

filter_order <- c(
  "NsLogit", "Logistic",
  "SCDE", "scDD", "MASTcpm",
  "MASTcpmDetRate", "ROTSvoom", "MASTtpm",
  "MASTtpmDetRate", "ROTScpm", "Wilcoxon",
  "ttest", "D3E", "limmatrend",
  "SeuratBimodIsExpr2", "metagenomeSeq", "edgeRQLFDetRate",
  "edgeRQLF", "DESeq2census", "ROTStpm",
  "voomlimma", "DEsingle", "SAMseq",
  "SeuratTobit", "edgeRLRTcensus", "monoclecensus",
  "edgeRLRTdeconv", "edgeRLRT", "monoclecount",
  "monocle", "DESeq2nofilt", "DESeq2",
  "NODES", "BPSC", "DESeq2betapFALSE",
  "edgeRLRTrobust", "SeuratBimod", "SeuratBimodnofilt")

res <- readRDS(paste0(work_dir, "sim123/sim123_res.rds"))
res <- rename(res, "LR", "Logistic") # rename LR to Logistic
method <- unique(res$FDR$Method)

# set conquer's method order
res <- order_method(res, unfilter_order)

nFac <- length(unique(method))
alpha <- rep(0, nFac)
alpha[1] <- 0.5
rec <- data.frame(xmin = head(seq <- seq(0.5, nFac + .5, 1), -1), xmax = tail(seq, -1), alpha = alpha)

fdr_png1 <- plot_all(res$FDR, "(a) FDR (unfiltered)", rec, NA) +
  geom_vline(xintercept = 0.05, color = "red", size = 0.6) +
  theme(legend.position = "none")
tpr_png1 <- plot_all(res$TPR, "(b) TPR (unfiltered)", rec, NA) +
  theme(legend.position = "none",
    axis.title.y=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank())
f1_png1 <- plot_all(res$F1, "(c) F1 (unfiltered)", rec, NA) +
  theme(legend.position = "none",
    axis.title.y=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank())

test <- plot_all(res$FDR, "", rec, NA) +
  theme(legend.position = "bottom",
    legend.box = "vertical")
legend <- cowplot::get_legend(test)

png(filename = paste0(work_dir, "Figure/no_filt.png"), width = 15, height = 15, units = 'in', res = 300)
grid.arrange(arrangeGrob(fdr_png1, tpr_png1, f1_png1, ncol = 3, widths=c(2, 1.2, 1.2)), legend, ncol=1, heights=c(10,1))
dev.off()

res <- readRDS(paste0(work_dir, "sim123_filtered/sim123_res_TPM_1_25p.rds"))
res <- rename(res, "LR", "Logistic") # rename LR to Logistic

# set conquer's method order
res <- order_method(res, filter_order)

# set highlight area
nFac <- length(unique(method))
alpha <- rep(0, nFac)
alpha[1] <- 0.5
rec <- data.frame(xmin = head(seq <- seq(0.5, nFac + .5, 1), -1), xmax = tail(seq, -1), alpha = alpha)

fdr_png2 <- plot_all(res$FDR, "(a) FDR (filtered)", rec, NA) +
  geom_vline(xintercept = 0.05, color = "red", size = 0.6) +
  theme(legend.position = "none")
tpr_png2 <- plot_all(res$TPR, "(b) TPR (filtered)", rec, NA) +
  theme(legend.position = "none",
    axis.title.y=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank())
f1_png2 <- plot_all(res$F1, "(c) F1 (filtered)", rec, NA) +
  theme(legend.position = "none",
    axis.title.y=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank())

png(filename = paste0(work_dir, "Figure/filt.png"), width = 15, height = 15, units = 'in', res = 300)
grid.arrange(arrangeGrob(fdr_png2, tpr_png2, f1_png2, ncol = 3, widths=c(2, 1.2, 1.2)), legend, ncol=1, heights=c(10,1))
dev.off()
