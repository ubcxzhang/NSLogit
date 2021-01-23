  
args = commandArgs(trailingOnly=TRUE)
work_dir <- as.character(args[1]) # the dir containing all files should not end with "/"
work_dir <- paste0(work_dir, "/")

# set R library
.libPaths(paste0(work_dir, "3.5"))

suppressPackageStartupMessages(library(ROTS))
suppressPackageStartupMessages(library(monocle))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(scde))
suppressPackageStartupMessages(library(DEsingle))
suppressPackageStartupMessages(library(BPSC))
suppressPackageStartupMessages(library(MAST))
suppressPackageStartupMessages(library(samr))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(metagenomeSeq))
suppressPackageStartupMessages(library(scDD))

sessionInfo()

