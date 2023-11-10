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
    opt_str = c("-o", "--output_sce_file"),
    type = "character",
    help = "Output path for clustered SCE file. Must end in .rds"
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

if (!stringr::str_ends(opt$output_sce_file, ".rds")) {
  stop("`output_sce_file` must end in .rds")
}

# only perform clustering if reduced dimension embeddings are present
# otherwise just return the object
if (!opt$pca_name %in% reducedDimNames(sce)) {
  warning("No reduced dimensions present with provided `pca_name`, skipping clustering")
} else {
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
  sce$cluster <- clusters
  # capitalize algorithm and weighting before adding to metadata
  metadata(sce)$cluster_algorithm <- stringr::str_to_sentence(opt$cluster_algorithm)
  metadata(sce)$cluster_weighting <- stringr::str_to_sentence(opt$cluster_weighting)
  metadata(sce)$cluster_nn <- opt$nearest_neighbors
}

# export -------------------
readr::write_rds(sce, opt$output_sce_file)
