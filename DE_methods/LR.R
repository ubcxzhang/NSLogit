suppressPackageStartupMessages(library(lmtest))

logit <- function(y, x) {
  table <- data.frame(x, y)
  full <-  glm(y~., data=table, family=binomial)
  null <-  glm(y~1, data=table, family=binomial)
  lrtest <- lrtest(full, null)
  lrtest$Pr[2]
}

p.logit <- function(y, xx) {
  # work_dir <- "/home/xsslnc/projects/def-ubcxzh/xsslnc/"
  # cl <- makeCluster(1)
  # clusterEvalQ(cl, {
  #  library(lmtest)
  #  source(work_dir, 'MasterProject/DE_methods/LR.R') # modify this
  # })
  # res <- parApply(cl, X=xx, MARGIN=2, FUN=logit, y=factor(y))
  # stopCluster(cl)
  res <- apply(X=xx, MARGIN=2, FUN=logit, y=factor(y))
  return(as.numeric(res))
}

run_LR <- function(data) {
  xx <- as.data.frame(t(data$count))
  y <- factor(data$condt)
  t <- system.time({
    res <- p.logit(y, xx)
  })

  list(
    time = t[3],
    df = data.frame(
      pval = res,
      padj = p.adjust(res, method="fdr", n=nrow(data$count)),
      row.names = rownames(data$count)))
}
