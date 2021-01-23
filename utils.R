suppressPackageStartupMessages(library(MultiAssayExperiment))
suppressPackageStartupMessages(library(survey))

#' Read rds file with name and path specified for dataset (for simulate data)
#' @param datapath string address of rds file
#' @param dataname stirng name of file without extension
#' @param name string name of dataset
#' @return a list of count matrix, tpm, group and dataset name
read_sim_data_from_rds <- function(datapath, dataname) {
  cat("[read_sim_data_from_rds] ")
  data <- readRDS(paste0(datapath,dataname,".rds"))
  data@sampleMap$assay <- factor(data@sampleMap$assay)
  data <- BiocGenerics::updateObject(data)
  #load(paste0(datapath,dataname,".rda"))
  count <- assays(experiments(data)[["gene"]])[["count_lstpm"]]
  tpm <- assays(experiments(data)[["gene"]])[["TPM"]]
  cat("Loading ", dataname, " ")
  cat(dim(count)[1], "genes ", dim(count)[2], "samples \n")

  # remove all 0 counts genes
  count <- count[rowSums(count) > 0, ]
  tpm <- tpm[rowSums(tpm) > 0, ]
  cat("After removing all 0 count genes ")
  cat(dim(count)[1], "genes ", dim(count)[2], "samples \n")

  test <- identical(rownames(colData(data)), colnames(count))&&identical(rownames(colData(data)), colnames(tpm))
  cat("colData(data)'s row matches count/tpm's column: ", paste0(test), "\n")
  condt <- colData(data)$group
  cat("Number of unique group: ", paste0(length(unique(condt))), "\n")
  cat("Considering groups:", paste(levels(factor(condt)), collapse = " vs "), "\n")

  list(count = count, tpm = tpm, condt = condt, name = dataname)
}

#' Remove ERCC spike-ins
#' @param data output of readRDS()
remove_ERCC <- function(data) {
  data <- subsetByRow(data, grep("^ERCC-", unique(unlist(rownames(data))),
    invert = TRUE, value = TRUE))
  data
}

#' Take one subset of data
take_sub <- function(data, subsets, sz, i) {
  keep_samples <- subsets$keep_samples
  imposed_condition <- subsets$out_condition
  s <- keep_samples[[as.character(sz)]][as.integer(i), ]
  count <- data$count[,s]
  tpm <- data$tpm[,s]
  condt <- imposed_condition[[as.character(sz)]][as.integer(i), ]
  names(condt) <- colnames(data$count[,s])
  sub_data <- list(count=count, tpm=tpm, condt=condt)
  sub_data
}

#' Read rds file for null data with name and path specified for dataset
#' @param datapath string address of rds file
#' @param dataname stirng name of file without extension
#' @param name string name of dataset
#' @param keepgroups string name of group to keep
#' @param clean bool whether need to remove ERCC spike-ins
#' @return a list of count matrix, tpm, group and dataset name
read_null_data_from_rds <- function(datapath, dataname, groupid, keepgroups, clean = FALSE) {
  data <- readRDS(paste0(datapath,dataname,".rds"))
  data@sampleMap$assay <- factor(data@sampleMap$assay)
  data <- BiocGenerics::updateObject(data)

  if (clean == TRUE) data <- remove_ERCC(data)

  count <- assays(experiments(data)[["gene"]])[["count_lstpm"]]
  tpm <- assays(experiments(data)[["gene"]])[["TPM"]]
  cat("Loading ", dataname, " ")
  cat(dim(count)[1], "genes ", dim(count)[2], "samples \n")

  # only keep one group
  test <- identical(rownames(colData(data)), colnames(count))&&identical(rownames(colData(data)), colnames(tpm))
  cat("Check colData(data)'s row matches count/tpm's column: ", paste0(test), "\n")

  if (length(groupid) > 1) {
    colData(data)[, paste(groupid, collapse = ".")] <-
      as.character(interaction(as.data.frame(colData(data)[, groupid])))
    groupid <- paste(groupid, collapse = ".")
  }

  keep <- colData(data)[,groupid] == keepgroups
  count <- count[, keep]
  tpm <- tpm[, keep]
  condt <- rep(keepgroups, dim(count)[2])
  condt <- factor(condt)
  cat("Keep group", keepgroups, "and use ")
  cat(dim(count)[1], "genes ", dim(count)[2], "samples \n")

  # remove all 0 counts genes in this group
  count <- count[rowSums(count) > 0, ]
  tpm <- tpm[rowSums(tpm) > 0, ]
  cat("After removing all 0 count genes ")
  cat(dim(count)[1], "genes ", dim(count)[2], "samples \n")

  list(count = count, tpm = tpm, condt = condt, name = dataname)
}

#' Generate subfile configuration
#' @param data list output from read_data_from_rds
#' @param seed int seed for sampling
#' @param sizes vector of size
#' @param nreps vector of numbers of subset for corresponding size
#' @param store_dir stirng address to store the ouput rds
#' @param filename stirng filename of output without extension
#'
gen_subfile <- function(data, seed, sizes, nreps, store_dir, filename) {
  # generate subfile for each dataset
  # 16 sub_datasets for each
  cat("[gen_subfile]: Generate subfile from size list", sizes, "and replicate list", nreps, "with seed", seed, "\n")
  condt <- data$condt
  names(condt) <- colnames(data$count) # cell's name
  names(sizes) <- sizes
  names(nreps) <- sizes
  set.seed(seed)

  if (length(unique(condt)) == 1) {
    condt <- condt[sort(sample(1:length(condt), 2 * max(sizes)))]
  } else {
    condt <- condt[sort(
      stratsample(
        as.character(condt),
        structure(rep(max(sizes), 2), names = levels(factor(condt)))
      )
    )]
  }

  ngroups <- nlevels(factor(condt))

  keep_tmp <- lapply(sizes, function(sz) {
    unique(t(sapply(1:nreps[as.character(sz)], function(i) {
      if (length(unique(condt)) == 1) {
        tmpn <- names(condt)
        condt2 <- paste0(condt, ".", sample(rep(c("1", "2"), ceiling(length(condt)/2)))[1:length(condt)])
        names(condt2) <- tmpn
        ngroups <- 2
      } else {
        condt2 <- condt
      }

      smp <- names(condt2)[sort(stratsample(as.character(condt2),
        structure(rep(sz, ngroups),
          names = levels(factor(condt2)))))]

      cdt <- condt2[smp]
      paste(smp, cdt, sep = "___")
    })))
  })

  keep_samples <- lapply(keep_tmp, function(w) {
    rbind(apply(w, 2, function(s) sapply(strsplit(s, "___"), .subset, 1)))})

  out_condition <- lapply(keep_tmp, function(w) {
    rbind(apply(w, 2, function(s) sapply(strsplit(s, "___"), .subset, 2)))})

  cat("The subfile is saved at",  paste0(store_dir, filename, ".rds \n"))

  saveRDS(list(keep_samples = keep_samples, out_condition = out_condition),
    file = paste0(store_dir, filename, ".rds"))
}

#' Apply DE approach on 16 sub dataset (with different size) of a simulated dataset
#' And write result in a rds
#' @param data list output from read_from_rds
#' @param use_method string the name of approach e.g. DESeq2
#' @param subfile string output from gen_subfile (name of subfile)
#' @param datapath string address of subfile and output file
#'
run_DE <- function(data, use_method, subfile, datapath) {
  subsets <- readRDS(paste0(datapath, subfile, ".rds"))
  keep_samples <- subsets$keep_samples
  imposed_condition <- subsets$out_condition
  sizes <- names(keep_samples)
  res <- list()
  cat("[run_DE]: Apply", use_method, "on", data$name, "... \n")

  for (sz in sizes) {
    cat("size:", sz, " ")
    # each rep of certain size
    n <- nrow(keep_samples[[as.character(sz)]])
    for (i in 1:n) {
      cat(i, " ")
      if (i == n) cat('\n')
      s <- keep_samples[[as.character(sz)]][i, ]
      count <- data$count[,s]
      tpm <- data$tpm[,s]
      condt <- imposed_condition[[as.character(sz)]][i, ]
      names(condt) <- colnames(data$count[,s])
      sub_data <- list(count=count, tpm=tpm, condt=condt)
      res[[paste0(use_method,'.',sz,'.',i)]] <- get(paste0("run_", use_method))(sub_data)
    }
  }

  # names(res) used_method.size.n_rep
  cat("Result saved at", paste0(datapath, data$name, "_", use_method, ".rds"))
  saveRDS(res, file = paste0(datapath, data$name, "_", use_method, ".rds"))
}

#' Apply DE approach on ny simulated dataset and write result in a rds
#' @param data string name of simulated dataset (without .rds)
#' @param use_method string the name of approach e.g. DESeq2
#' @param datapath string address of dataset file
#' @param storepath string address of output file
#'
run_DE_sim <- function(data, use_method, datapath, savepath) {
  sim <- readRDS(paste0(datapath, data, ".rds"))
  res <- list()
  res[[paste0(use_method)]] <- get(paste0("run_", use_method))(sim)
  saveRDS(res, file = paste0(savepath, data, "_", use_method, ".rds"))
}
