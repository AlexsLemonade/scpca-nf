#!/usr/bin/env Rscript

# Script used to perform clustering on a given SCE object
# 
# This script reads in an RDS file containing an SCE, and performed
#  graph-based clustering using the specified algorithm. Cluster identities are 
#  stored in the SCE's colData slot, and the SCE is written back out to the
#  original RDS file.

# import libraries
library(optparse)
library(SingleCellExperiment)

# Set up optparse options
option_list <- list(
  make_option(
    opt_str = c("-i", "--sce_file"),
    type = "character",
    help = "Path to RDS file that contains the SCE object to cluster."
  ),
  make_option(
    opt_str = c("--pca_name"),
    type = "character",
    default = "PCA",
    help = "Name of the PCA reduced dimensionality representation to perform
      clustering on."
  ), 
  make_option(
    opt_str = c("--cluster_algorithm"),
    type = "character",
    default = "leiden",
    help = "Clustering algorithm to use. Must be one of the options available in bluster.
      Default is 'leiden'."
  ),
  make_option(
    opt_str = c("--nn"),
    type = "integer",
    default = 15, 
    help = "Nearest neighbors parameter to set for graph-based clustering."
  ),
  make_option(
    opt_str = c("--seed"),
    type = "integer",
    help = "random seed to set during clustering."
  )
)

# Setup ------------------------------
# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

set.seed(opt$seed)

# check and read in SCE file
if (!file.exists(opt$sce_file)) {
  stop("Input `sce_file` is missing.")
}
sce <- readr::read_rds(opt$sce_file)


# check pca_name is present
if (!opt$pca_name %in% reducedDimNames(sce)) {
  stop("Provided `pca_name` is not present in the SCE object.")
}

# Perform clustering ----------------

# extract the principal components matrix
pca_matrix <- reducedDim(sce, opt$pca_name)

# cluster with specified algorithm
clusters <- bluster::clusterRows(
  pca_matrix,
  bluster::SNNGraphParam(
    k = opt$nn,
    cluster.fun = opt$cluster_algorithm
  )
)

# add clusters and associated parameters to SCE object
sce$clusters <- clusters
metadata(sce)$cluster_algorithm <- opt$cluster_algorithm
metadata(sce)$cluster_nn <- opt$nn


# export -------------------
readr::write_rds(sce, opt$sce_file)
