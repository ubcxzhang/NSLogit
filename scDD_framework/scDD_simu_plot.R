args = commandArgs(trailingOnly=TRUE)
work_dir <- as.character(args[1]) # the dir containing all files should not end with "/"
work_dir <- paste0(work_dir, "/")

# set R library
.libPaths(paste0(work_dir, "3.5"))

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(ggplot2))

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


# load file
scDD_simu_res <- readRDS(paste0(work_dir, "scDD_simu_result/comb_results/scDD_simu_res.rds"))
FDR <- scDD_simu_res$FDR
TPR <- scDD_simu_res$TPR
DD  <- scDD_simu_res$DD
F1  <- scDD_simu_res$F1

# rename logistic
FDR[FDR$Method == "LR", "Method"] <- "Logistic"
TPR[TPR$Method == "LR", "Method"] <- "Logistic"
F1[F1$Method == "LR", "Method"] <- "Logistic"
DD[DD$Method == "LR", "Method"] <- "Logistic"

# generate shaded area information
nFac <- length(unique(FDR$Method))
alpha <- rep(0, nFac)
alpha[19] <- 0.5
rec <- data.frame(xmin = head(seq <- seq(0.5, nFac + .5, 1), -1), xmax = tail(seq, -1), alpha = alpha)

# get plot
fdr_png <- plot_all(FDR, "(a) FDR", rec, NA) +
  geom_vline(xintercept = 0.05, color = "red", size = 0.6) +
  theme(legend.position = "none")
tpr_png <- plot_all(TPR, "(b) TPR", rec, NA) +
  theme(legend.position = "none",
    axis.title.y=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank())
f1_png <- plot_all(F1, "(c) F1", rec, NA) +
  theme(legend.position = "none",
    axis.title.y=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank())

de_png <- plot_all(DD %>% filter(Gene=="DE"), "(a) DE", rec, NA) +
  theme(legend.position = "none")
db_png <- plot_all(DD %>% filter(Gene=="DB"), "(d) DB", rec, NA) +
  theme(legend.position = "none",
    axis.title.y=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank())
dm_png <- plot_all(DD %>% filter(Gene=="DM"), "(c) DM", rec, NA) +
  theme(legend.position = "none",
    axis.title.y=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank())
dp_png <- plot_all(DD %>% filter(Gene=="DP"), "(b) DP", rec, NA) +
  theme(legend.position = "none",
    axis.title.y=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank())

test <- plot_all(TPR, "(b) TPR", rec, NA) +
  theme(legend.position = "bottom",
        legend.box = "vertical")
legend <- cowplot::get_legend(test)

png(filename = paste0(work_dir, "Figure/scDD_result2.png"), width = 15, height = 15, units = 'in', res = 300)
grid.arrange(arrangeGrob(de_png, dp_png, dm_png, db_png, ncol = 4, widths=c(2, 1.1, 1.1, 1.1)), legend, ncol=1, heights=c(10,1))
dev.off()

png(filename = paste0(work_dir, "Figure/scDD_result1.png"), width = 15, height = 15, units = 'in', res = 300)
grid.arrange(arrangeGrob(fdr_png, tpr_png, f1_png, ncol = 3, widths=c(2, 1.2, 1.2)), legend, ncol=1, heights=c(10,1))
dev.off()

