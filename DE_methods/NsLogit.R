suppressPackageStartupMessages(library(splines))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(glmnet))

fit_bic <- function(fit) {
  tLL <- deviance(fit)
  k <- fit$df + 1
  n <- fit$nobs
  return(log(n)*k + tLL)
}

lr_p <- function(fit) {
  df <- fit$df
  if (df == 0) {return(1)}
  test <- fit$nulldev - deviance(fit)
  return(pchisq(test, df=df, lower.tail=FALSE))
}

# Spline Knots
inner.knots<- function(xx, K)
{
  uxx <- unique(xx)
  if (length(uxx) < 2) return(rep(uxx, K))
  tmp1 <- weighted.mean(head(sort(uxx),2),c(0.1,0.9))
  if (K == 1) return(tmp1)
  tmp2 <- xx[xx != min(xx)]
  return(c(tmp1, quantile(tmp2, probs = 1:(K-1)/K)))
}

separation_chi_test <- function(y, x, e) {
  grp1 <- x[y==unique(y)[1]]; grp2 <- x[y==unique(y)[2]]
  if (sum(grp1)!=0&sum(grp2)!=0) {
    return(NA) # no perfect or quasi-separation
  } else {
    # number of sample in each group
    ng1 = length(grp1); ng2 = length(grp2)

    g1_0 = sum(as.numeric(grp1==0)) + e
    g2_0 = sum(as.numeric(grp2==0)) + e
    g1_1 = sum(as.numeric(grp1!=0)) + e
    g2_1 = sum(as.numeric(grp2!=0)) + e

    e1_0 = (g1_0 + g2_0)*ng1/(ng1 + ng2)
    e1_1 = (g1_1 + g2_1)*ng1/(ng1 + ng2)

    e2_0 = (g1_0 + g2_0)*ng2/(ng1 + ng2)
    e2_1 = (g1_1 + g2_1)*ng2/(ng1 + ng2)

    test = (g1_0-e1_0)^2/e1_0 + (g2_0-e2_0)^2/e2_0 +
      (g1_1-e1_1)^2/e1_1 + (g2_1-e2_1)^2/e2_1
    return(pchisq(test, df=1, lower.tail = FALSE))
  }
}

ns.logit <- function(y, x, df=3)
{
  chi_p <- separation_chi_test(y, x, 0)
  if (is.na(chi_p)==FALSE) {return(chi_p)}

  if(length(x)/5 < df) {
    df <- floor(length(x)/5)
    cat(paste("df is too big, changed to", floor(length(x)/5)))
  }

  if (df>1){
    models <- vector("list",df-1)
    BICs <- rep(Inf, df-1)

    for (ii in 2:df) {
      if (length(unique(x))>ii) myknots <- inner.knots(x, K=ii-1) else myknots <- NULL
      dat <- try(data.frame(splines::ns(x, df=ii, knots=myknots)), silent=T)
      if(class(dat)=="try-error") break
      fit <- glmnet(as.matrix(dat), y, family="binomial", alpha=1, lambda=1/nrow(dat))
      models[[ii-1]] <- fit
      BICs[ii-1] <- fit_bic(fit)
    }
    # model comparisions
    bestM <- which.min(BICs)
    if (BICs[bestM]!=Inf) {
      pval <- lr_p(models[[bestM]])
    } else {
      pval <- NA
    }
  } else {
    model0 <- glm(y ~ 1, family="binomial")
    model  <- glm(y ~ x, family="binomial")
    pval <- anova(model, model0, test='Chisq')$`Pr(>Chi)`[2]
  }

  return(pval)
}

p.ns.logit <- function(y, xx, df=3) {
  #cl <- makeCluster(1)
  #clusterEvalQ(cl, {
  #  library(splines)
  #  source('/home/xsslnc/Spline/conquer_script/NsLogit.R') # modify this
  #})
  #res <- parApply(cl, X=xx, MARGIN=2, FUN=ns.logit, y=factor(y), df=df)
  #stopCluster(cl)
  res <- apply(X=xx, MARGIN=2, FUN=ns.logit, y=factor(y), df=df)
  return(as.numeric(res))
}

run_NsLogit <- function(data) {
  xx <- as.data.frame(t(data$count))
  y <- factor(data$condt)
  res <- p.ns.logit(y, xx, df=3)
  list(
    df = data.frame(
      pval = res,
      padj = p.adjust(res, method="fdr", n=nrow(data$count)),
      row.names = rownames(data$count)))
}
