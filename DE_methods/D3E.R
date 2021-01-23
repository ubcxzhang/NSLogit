run_D3E <- function(data) {
  tmp <- cbind(GeneID = rownames(data$count), data$count)
  colnames(tmp) <- c("GeneID", make.names(data$condt))
  rnb <- paste0(format(Sys.time(), "%Y_%m_%d_%H_%M_%S"), "_", round(runif(1) * 1e8))
  write.table(tmp, file = paste0(work_dir, "tmp/d3e_", rnb, ".txt"), row.names = FALSE,
    col.names = TRUE, quote = FALSE, sep = "\t")
  setwd(work_dir)
  cmd <- sprintf("python D3E/D3ECmd.py %s %s %s %s -m 0 -t 0 -z 0 -n 1 -v",
    paste0("tmp/d3e_", rnb, ".txt"),
    paste0("tmp/d3e_", rnb, ".out"),
    levels(factor(make.names(data$condt)))[1],
    levels(factor(make.names(data$condt)))[2])
  system(cmd)
  res <- read.delim(paste0(work_dir, "tmp/d3e_", rnb, ".out"), header = TRUE,
    as.is = TRUE, row.names = 1)

  list(df = data.frame(
    pval = res$p.value,
    padj = p.adjust(res$p.value, method = "BH"),
    row.names = rownames(res)))
}


