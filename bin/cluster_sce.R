#!/usr/bin/env Rscript

# Script used to perform clustering on a given SCE object
# 
# This script reads in an RDS file containing an SCE, and performs
#  graph-based clustering using the specified algorithm. Cluster identities are 
#  stored in the SCE's colData slot, and the SCE is written back out to the
#  original RDS file.

# import libraries
library(optparse)
library(SingleCellExperiment)

# Set up optparse options
option_list <- list(
  make_option(
    opt_str = c("-i", "--input_sce_file"),
    type = "character",
    help = "Path to RDS file that contains the processed sce SCE object to cluster."
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
    default = "louvain",
    help = "Clustering algorithm to use. Must be one of the options available in bluster.
      Default is 'louvain'."
  ),
  make_option(
    opt_str = c("--nearest_neighbors"),
    type = "integer",
    default = 15, 
    help = "Nearest neighbors parameter to set for graph-based clustering."
  ),
  make_option(
    opt_str = c("--random_seed"),
    type = "integer",
    help = "random seed to set during clustering."
  )
)

# Setup ------------------------------
# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

set.seed(opt$seed)

# check and read in SCE file
if (!file.exists(opt$input_sce_file)) {
  stop("Input `input_sce_file` is missing.")
}
sce <- readr::read_rds(opt$input_sce_file)


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
  bluster::NNGraphParam(
    k = opt$nearest_neighbors,
    type = "jaccard",
    cluster.fun = opt$cluster_algorithm
  )
) 
# TODO: may wish to modify additional cluster parameters depending on cluster algorithm

# add clusters and associated parameters to SCE object
sce$clusters <- clusters
metadata(sce)$cluster_algorithm <- opt$cluster_algorithm
metadata(sce)$cluster_nn <- opt$nearest_neighbors

# export -------------------
# we are overwriting the input file here:
readr::write_rds(sce, opt$input_sce_file)
