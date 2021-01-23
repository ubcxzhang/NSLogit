avglogfc_pct <- function(
  data,
  assay,
  slot,
  thresh.min = 0,
  min.pct = 0.1,
  min.diff.pct = -Inf,
  logfc.threshold = 0.25,
  min.cells.feature = 3) {

  library(Seurat)
  library(dplyr)
  library(MASS)

  Idents(object = data) <- "Status"
  data.use <- GetAssayData(object = data[[assay]], slot = slot)
  data.use <- as.matrix(data.use)

  cells.1 <- WhichCells(data, idents = "ILD")
  cells.2 <- WhichCells(data, idents = "Control")

  features <- rownames(data.use)

  pct.1 <- round(
    x = rowSums(x = data.use[features, cells.1, drop = FALSE] > thresh.min) /
      length(x = cells.1),
    digits = 3
  )
  pct.2 <- round(
    x = rowSums(x = data.use[features, cells.2, drop = FALSE] > thresh.min) /
      length(x = cells.2),
    digits = 3
  )

  data.alpha <- cbind(pct.1, pct.2)
  colnames(x = data.alpha) <- c("pct.1", "pct.2")

  # filter out gene if both pcts are less than min.pct
  alpha.min <- apply(X = data.alpha, MARGIN = 1, FUN = max)
  names(x = alpha.min) <- rownames(x = data.alpha)
  features <- names(x = which(x = alpha.min > min.pct))
  alpha.diff <- alpha.min - apply(X = data.alpha, MARGIN = 1, FUN = min)

  # filter out gene if absolute difference of pct is less than min.diff.pct
  features <- names(
    x = which(x = alpha.min > min.pct & alpha.diff > min.diff.pct)
  )

  pseudocount.use <- 1

  data.1 <- apply(
    X = data.use[features, cells.1, drop = FALSE],
    MARGIN = 1,
    FUN = function(x) {
      return(log(x = mean(x = expm1(x = x)) + pseudocount.use))
    }
  )

  data.2 <- apply(
    X = data.use[features, cells.2, drop = FALSE],
    MARGIN = 1,
    FUN = function(x) {
      return(log(x = mean(x = expm1(x = x)) + pseudocount.use))
    }
  )
  total.diff <- (data.1 - data.2)

  # filter out gene if absolute log fold change is less than logfc.threshold
  features.diff <- names(x = which(x = abs(x = total.diff) > logfc.threshold))
  features <- intersect(x = features, y = features.diff)

  cat(paste0("selected genes for DE test: ", length(features), "\n"))
  select_ngene <- features

  ngene_before_filter <- nrow(data.use)

  data.use = data.use[features, c(cells.1, cells.2), drop = FALSE]

  de.results <- GLMDETest(data.use, cells.1, cells.2, min.cells = min.cells.feature)
  de.results[, "avg_logFC"] <- total.diff[rownames(x = de.results)]
  de.results <- cbind(de.results, data.alpha[rownames(x = de.results), , drop = FALSE])

  de.results$p_val_adj = p.adjust(
    p = de.results$p_val,method = "bonferroni",
    n = ngene_before_filter
  )
  return(de.results)
}

GLMDETest <- function(
  data.use,
  cells.1,
  cells.2,
  min.cells = 3,
  latent.vars = NULL
) {
  group.info <- data.frame(
    group = rep(
      x = c('Group1', 'Group2'),
      times = c(length(x = cells.1), length(x = cells.2))
    )
  )
  rownames(group.info) <- c(cells.1, cells.2)
  group.info[, "group"] <- factor(x = group.info[, "group"])
  latent.vars <- if (is.null(x = latent.vars)) {
    group.info
  } else {
    cbind(x = group.info, latent.vars)
  }
  latent.var.names <- colnames(x = latent.vars)

  p_val <- rep(NA, nrow(data.use))
  ncell_gp1 <- rep(NA, nrow(data.use))
  ncell_gp2 <- rep(NA, nrow(data.use))
  var_gp <- rep(NA, nrow(data.use))

  for (i in 1:nrow(data.use)) {
    latent.vars[, "GENE"] <- as.numeric(x = data.use[i, ])
    ncell_gp1[i] <- sum(latent.vars$GENE[latent.vars$group == "Group1"] > 0)
    ncell_gp2[i] <- sum(latent.vars$GENE[latent.vars$group == "Group2"] > 0)
    var_gp[i] <- var(x = latent.vars$GENE)

    fmla <- as.formula(object = paste(
      "GENE ~",
      paste(latent.var.names, collapse = "+")
    ))

    p.estimate <- NA
    try(
      expr = p.estimate <- summary(
        object = glm.nb(formula = fmla, data = latent.vars)
      )$coef[2, 4],
      silent = TRUE
    )
    p_val[i] <- p.estimate
  }

  features.keep <- rownames(data.use)
  to.return <- cbind(p_val, ncell_gp1, ncell_gp2, var_gp)
  rownames(to.return) <- features.keep
  colnames(to.return) <- c("p_val", "ncells_in_group1", "ncells_in_group2", "variance")
  return(as.data.frame(to.return))
}
