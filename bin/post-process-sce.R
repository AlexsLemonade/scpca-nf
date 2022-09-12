#!/usr/bin/env Rscript

# This script takes a SingleCellExperiment stored in a .rds file with empty droplets removed 
# and removes low quality cells, performs normalization, and dimensionality reduction.

# import libraries
suppressPackageStartupMessages({
  library(optparse)
  library(SingleCellExperiment)
})
# set up arguments
option_list <- list(
  make_option(
    opt_str = c("--library_id"),
    type = "character",
    help = "library id"
  ),
  make_option(
    opt_str = c("-u", "--input_sce_file"),
    type = "character",
    help = "path to rds file with input sce object to be processed"
  ),
  make_option(
    opt_str = c("-f", "--output_sce_file"),
    type = "character",
    help = "path to output rds file to store processed sce object. Must end in .rds"
  ),
  make_option(
    opt_str = c("--prob_compromised_cutoff"),
    type = "double",
    default = 0.75,
    help = "probability compromised cutoff used for filtering cells with miQC",
    metavar = "double"
  ),
  make_option(
    opt_str = c("--gene_cutoff"),
    type = "integer",
    default = 200,
    help = "Minimum number of genes per cell cutoff used for filtering cells."
  ),
  make_option(
    opt_str = c("-n", "--top_n"),
    type = "double",
    default = 2000,
    help = "top number of high variance genes to use for dimension reduction;
            the default is top_n = 2000",
  ),
  make_option(
    opt_str = c("-r", "--random_seed"),
    type = "integer",
    help = "A random seed for reproducibility."
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

# set seed
set.seed(opt$random_seed)

# check that input file file exists
if(!file.exists(opt$input_sce_file)){
  stop("Missing unfiltered.rds file")
}

# check that output file name ends in .rds
if(!(stringr::str_ends(opt$output_sce_file, ".rds"))){
  stop("output file name must end in .rds")
}

# read in input rds file
sce <- readr::read_rds(opt$input_sce_file)

# check that prob compromised cutoff is between 0-1
if(!dplyr::between(opt$prob_compromised_cutoff, 0, 1)){
  stop("--prob_compromised_cutoff must be a number between 0 to 1")
}

# create ccdl_filter column
if(all(is.na(sce$prob_compromised))){
  colData(sce)$ccdl_filter <- ifelse(
    sce$detected >= opt$gene_cutoff, "Keep", "Remove"
  )
  metadata(sce)$ccdl_filter_method <- "Minimum_gene_cutoff"
} else {
  # remove cells with >= probability compromised cutoff + min gene cutoff
  colData(sce)$ccdl_filter <- ifelse(
    sce$prob_compromised < opt$prob_compromised_cutoff & 
      sce$detected >= opt$gene_cutoff, 
    "Keep", 
    "Remove"
  )
  metadata(sce)$ccdl_filter_method <- "miQC"
}

# filter sce using criteria in ccdl_filter 
filtered_sce <- sce[, sce$ccdl_filter == "Keep"]

# cluster prior to normalization 
qclust <- NULL
tryCatch({
  # Cluster similar cells
  qclust <- scran::quickCluster(filtered_sce) 
  
  # Compute sum factors for each cell cluster grouping
  filtered_sce <- scran::computeSumFactors(filtered_sce, clusters = qclust)
  
  # Include note in metadata re: clustering before computing sum factors
  metadata(filtered_sce)$normalization <- "deconvolution"
})

if (is.null(qclust)) {
  
  # Include note in metadata re: failed clustering
  metadata(filtered_sce)$normalization <- "log-normalization"
  
  # Keep positive counts for `logNormCounts()`
  filtered_sce <- filtered_sce[, colSums(counts(filtered_sce)) > 0]
  
}

# Normalize and log transform
normalized_sce <- scater::logNormCounts(filtered_sce)

# remove filtered SCE to save space 
rm(filtered_sce)

# model gene variance using `scran:modelGeneVar()`
gene_variance <- scran::modelGeneVar(normalized_sce)

# select the most variable genes
var_genes <- scran::getTopHVGs(gene_variance, n = opt$top_n)

# save the most variable genes to the metadata
metadata(normalized_sce)$variable_genes <- var_genes

# dimensionality reduction 
# highly variable genes are used as input to PCA 
normalized_sce <- scater::runPCA(normalized_sce, 
                                 subset_row = var_genes)

# calculate a UMAP matrix using the PCA results
normalized_sce <- scater::runUMAP(normalized_sce, 
                                  dimred = "PCA")

# write out original SCE with additional filtering column
readr::write_rds(sce, opt$input_sce_file)

# write out processed SCE 
readr::write_rds(normalized_sce, opt$output_sce_file)

