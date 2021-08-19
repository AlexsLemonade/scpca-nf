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

# write filtered sce to output
readr::write_rds(filtered_sce, opt$filtered_file)

