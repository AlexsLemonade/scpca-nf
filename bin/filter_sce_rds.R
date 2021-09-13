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

# need to remove old gene-level rowData first 
rowData(filtered_sce) <- NULL

# recalculate rowData and add to filtered sce 
filtered_sce <- filtered_sce %>%
  scater::addPerFeatureQC()
  
# add prob_compromised to colData from miQC::mixtureModel 
model <- miQC::mixtureModel(filtered_sce)
filtered_sce <- miQC::filterCells(filtered_sce, model, posterior_cutoff = 1, verbose = FALSE)

# grab names of altExp, if any
alt_names <- altExpNames(filtered_sce)

for (alt in alt_names) {
  # remove old row data from unfiltered 
  rowData(altExp(test, alt)) <- NULL
  
  # add alt experiment features stats for filtered data
  altExp(test, alt) <- scater::addPerFeatureQC(altExp(test, alt))
}

# write filtered sce to output
readr::write_rds(filtered_sce, opt$filtered_file, compress = "gz")

