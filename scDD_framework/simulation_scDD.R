suppressPackageStartupMessages(library(MultiAssayExperiment))
suppressPackageStartupMessages(library(survey))
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(scDD))

#' format conquer dataset for scDD simulation input requirement
#' @param datapath string address of original dataset rds file
#' @param dataname string name of original dataset rds file
#' @param groupid the colname in colData(data) to use for group identification
#' @param keepgroup a list of group name (length = 2) to use for simulation
#' @return a singlecellexperiment object
#'
format_scDD <- function(datapath, dataname, groupid, keepgroup) {
  data <- readRDS(paste0(datapath, dataname, ".rds"))
  data@sampleMap$assay <- factor(data@sampleMap$assay)
  data <- BiocGenerics::updateObject(data)

  count <- assays(experiments(data)[["gene"]])[["count_lstpm"]]
  cat("Loading ", dataname, " ")
  cat(dim(count)[1], "genes ", dim(count)[2], "samples \n")

  # check cells order in colData (row) and count/tpm (column)
  test <- identical(rownames(colData(data)), colnames(count))
  cat("Check colData(data)'s row matches count's column: ", paste0(test), "\n")

  if (length(groupid) > 1) {
    colData(data)[, paste(groupid, collapse = ".")] <-
      as.character(interaction(as.data.frame(colData(data)[, groupid])))
    groupid <- paste(groupid, collapse = ".")
  }

  # keep cells of keepgroup
  keep <- colData(data)[, groupid] %in% keepgroup
  count <- count[, keep]
  condt <- colData(data)[, groupid][keep]
  cat("Considering groups:", paste(levels(factor(condt)), collapse = " vs "), "\n")
  cat("Using Cells: ", dim(count)[2], "\n")

  # condition should be 1 or 2
  condition <- condt
  condition <- as.character(condition)
  condition[condition==unique(condition)[1]] <- 1
  condition[condition==unique(condition)[2]] <- 2
  levels(condition)[1] <- 1
  levels(condition)[2] <- 2
  condition <- as.numeric(condition)
  names(condition) <- colnames(count)

  # create SingleCellExperiment Object
  sce <- SingleCellExperiment(
    assays = list(counts=count), # the normcounts assays slot contains one matrix
    colData = data.frame(condition))
  return(sce)
}

#' simulation using scDD: 90% of non-single gene
#' @param sce output of format_scDD
#' @param nsample integer number of cells per gene
#' @param seed simulation seed
#' @return a singlecellexperiment object
simulation_scDD <- function(sce, nsample, seed) {
  set.seed(12345)
  ngene <- dim(sce)[1] # simulate same number of genes as original dataset

  # equal number of four genes with signal
  nDE <- nDP <- nDM <- nDB <- round(ngene*0.1*0.25)
  nEE <- nEP <- round(ngene*0.9*0.5)

  # screen and normalization (genes with > 90% zeros are removed)
  sce.scran <- preprocess(sce, zero.thresh=0.9, scran_norm=TRUE)

  # simulation: non-parallel
  SD <- simulateSet(
    sce.scran, numSamples=nsample,
    nDE=nDE, nDP=nDP, nDM=nDM, nDB=nDB,
    nEE=nEE, nEP=nEP, plots=FALSE,
    random.seed=seed,
    param=BiocParallel::MulticoreParam(workers=1))

  cat("After removing genes with > 90% zeros :", dim(counts(sce.scran))[1], "genes ", dim(counts(sce.scran))[2], "cells for simulation. \n")
  cat("Data is normalized\n")
  cat(paste0("Simulating ", nsample, " Cells...\n"))
  cat(paste0("Using seed ", seed, "\n"))
  cat(paste0("Significant genes ", nDE, " per types (DE, DP, DM, DB)\n"))
  cat(paste0("Rest of genes ", nEE, " per types (EP, EE)\n"))
  cat("Done!\n")

  return(SD)
}
