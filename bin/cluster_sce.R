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
    opt_str = c("--processed_sce_file"),
    type = "character",
    help = "Path to RDS file that contains the processed SCE object to cluster.
      Must contain a PCA matrix to calculate clusters from."
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
    opt_str = c("--cluster_weighting"),
    type = "character",
    default = "jaccard",
    help = "The type of weighting scheme to use for shared neighbors when performing
      graph-based clustering. Default is 'jaccard'."
  ),
  make_option(
    opt_str = c("--nearest_neighbors"),
    type = "integer",
    default = 20,
    help = "Nearest neighbors parameter to set for graph-based clustering. Default is 20."
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
if (!file.exists(opt$processed_sce_file)) {
  stop("Input `processed_sce_file` is missing.")
}
sce <- readr::read_rds(opt$processed_sce_file)


# check pca_name is present
if (!opt$pca_name %in% reducedDimNames(sce)) {
  stop("Provided `pca_name` is not present in the SCE object.")
}

# Perform clustering ----------------

# extract the principal components matrix
clusters <- scran::clusterCells(
  sce,
  use.dimred = opt$pca_name,
  BLUSPARAM = bluster::NNGraphParam(
    k = opt$nearest_neighbors,
    type = opt$cluster_weighting,
    cluster.fun = opt$cluster_algorithm
  )
)

# add clusters and associated parameters to SCE object
sce$clusters <- clusters
metadata(sce)$cluster_algorithm <- opt$cluster_algorithm
metadata(sce)$cluster_weighting <- opt$cluster_weighting
metadata(sce)$cluster_nn <- opt$nearest_neighbors

# export -------------------
# we are overwriting the `processed_sce_file` here
readr::write_rds(sce, opt$processed_sce_file)
