#!/usr/bin/env Rscript

# This script takes the output folder from alevin-fry as input and
# returns the filtered counts matrices as a SingleCellExperiment stored in a .rds file

# import libraries
library(magrittr)
library(optparse)

# set up arguments
option_list <- list(
  make_option(
    opt_str = c("-u", "--unfiltered_file"),
    type = "character",
    help = "path to rds file with unfiltered sce object"
  ),
  make_option(
    opt_str = c("-f", "--filtered_file"),
    type = "character",
    help = "path to output filtered rds file. Must end in .rds"
  ),
  make_option(
    opt_str = c("-m", "--mito_file"),
    type = "character",
    default = "",
    help = "path to list of mitochondrial genes"
  ), 
  make_option(
    opt_str = c("-n", "--feature_name"),
    type = "character",
    default = "ALT",
    help = "Feature type"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

# check that unfiltered file file exists
if(!file.exists(opt$unfiltered_file)){
  stop("Missing unfiltered.rds file")
}

# check that output file name ends in .rds
if(!(stringr::str_ends(opt$filtered_file, ".rds"))){
  stop("filtered file name must end in .rds")
}

# read in unfiltered rds file
unfiltered_sce <- readr::read_rds(opt$unfiltered_file)

# filter sce
filtered_sce <- scpcaTools::filter_counts(unfiltered_sce)

# need to remove old colData and rowData first 
colData(filtered_sce) <- NULL
rowData(filtered_sce) <- NULL

# add colData and rowData to filtered sce 
filtered_sce <- filtered_sce %>%
  scpcaTools::add_cell_mito_qc(mito = mito_genes, miQC = TRUE) %>%
  scater::addPerFeatureQC()

# if altExp is present, add feature data
if (opt$feature_name != "") {
  # remove old row data from unfiltered 
  rowData(altExp(filtered_sce, opt$feature_name)) <- NULL
  
  # add alt experiment features stats for filtered data
  altExp(filtered_sce, opt$feature_name) <- scater::addPerFeatureQC(altExp(filtered_sce, opt$feature_name))
}

# write filtered sce to output
readr::write_rds(filtered_sce, opt$filtered_file, compress = "gz")

