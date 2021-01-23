suppressPackageStartupMessages(library(scde))

run_scde <- function(data) {

  intcount <- apply(round(data$count), 2, function(x) {storage.mode(x) <- 'integer'; x})

  o.ifm <- scde.error.models(
    counts = intcount, groups = data$condt, n.cores = 1,
    threshold.segmentation = TRUE,
    save.crossfit.plots = FALSE, save.model.plots = FALSE,
    verbose = 0, min.size.entries = min(2000, nrow(data$count) - 1))

  valid.cells <- o.ifm$corr.a > 0

  o.ifm <- o.ifm[valid.cells, ]
  o.prior <- scde.expression.prior(models = o.ifm, counts = intcount[, valid.cells],
    length.out = 400, show.plot = FALSE)

  grp <- factor(data$condt[which(valid.cells)])
  names(grp) <- rownames(o.ifm)
  ediff <- scde.expression.difference(o.ifm, intcount[, valid.cells], o.prior,
    groups = grp, n.randomizations = 100,
    n.cores = 1, verbose = 0)
  p.values <- 2*pnorm(abs(ediff$Z), lower.tail = FALSE)
  p.values.adj <- 2*pnorm(abs(ediff$cZ), lower.tail = FALSE)

  list(df =
    data.frame(
      pval = p.values,
      padj = p.values.adj,
      row.names = rownames(ediff)))
}