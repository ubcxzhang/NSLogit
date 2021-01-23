args = commandArgs(trailingOnly=TRUE)
work_dir <- as.character(args[1]) # the dir containing all files should not end with "/"
work_dir <- paste0(work_dir, "/")

# set R library
.libPaths(paste0(work_dir, "3.5"))

suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))

# Plot time
scDD_simu_res <- readRDS(paste0(work_dir, "scDD_simu_result/comb_results/scDD_simu_res.rds"))
FDR <- scDD_simu_res$FDR
FDR[FDR$Method == "LR", "Method"] <- "Logistic"

time <- data.frame()
for (mt in unique(FDR$Method)) {
  for (sz in unique(FDR$Size)) {
    select <- FDR%>%filter(Method==mt&Size==sz)
    if (nrow(select)==0) {next()}
    avg <- mean(select$Time)
    res <- data.frame()
    res[1, "Method"] <- mt; res[1, "Time"] <- avg; res[1, "Size"] <- sz
    time <- rbind(time, res)
  }
}

time$Size <- as.character(time$Size)
time$Size = factor(time$Size, levels=
    c("50", "100", "500", "1000", "1500", "2000", "2500", "3000", "3500")
)
time$Method = factor(time$Method, levels=
    c(unique(time$Method)[unique(time$Method) != 'NsLogit'], "NsLogit")
)
col <- rep("grey", length(unique(time$Method)))
names(col) <- unique(time$Method)
col["NsLogit"] <- "red"

p <- ggplot(time, aes(x = Size, y = Time, group = Method, color = Method)) +
  geom_line() +
  scale_color_manual(values = col) +
  labs(x = "number of cells per group", y = "average runtime (second)", size=25) +
  theme(
    legend.position = "none",
    panel.background = element_rect(fill = "white", colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())

png(filename = paste0(work_dir, "Figure/scDD_time.png"), width = 5, height = 5, units = 'in', res = 300)
grid.arrange(p, ncol = 1)
dev.off()

cat(paste0("scDD_time.png saved at ", work_dir, "Figure\n"))
