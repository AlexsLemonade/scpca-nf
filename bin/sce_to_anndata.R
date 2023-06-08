#!/usr/bin/env Rscript

# This script takes a SingleCellExperiment stored in a .rds file and converts it
# to an AnnData object saved as an hdf5 file

# import libraries
suppressPackageStartupMessages({
  library(optparse)
  library(SingleCellExperiment)
})
# set up arguments
option_list <- list(
  make_option(
    opt_str = c("-i", "--input_sce_file"),
    type = "character",
    help = "path to rds file with input sce object to be converted"
  ),
  make_option(
    opt_str = c("-o", "--output_h5_file"),
    type = "character",
    help = "path to output hdf5 file to store AnnData object. Must end in .hdf5 or .h5"
  )
)

# Parse and apply input options -------------
opt <- parse_args(OptionParser(option_list = option_list))

# check that filtered SCE file exists
if(!file.exists(opt$input_sce_file)){
  stop(glue::glue("{opt$input_sce_file} does not exist."))
}

# check that output file name ends in .rds
if(!(stringr::str_ends(opt$output_h5_file, ".hdf5|.h5"))){
  stop("output file name must end in .hdf5 or .h5")
}

sce <- readr::read_rds(opt$input_sce_file)
scpcaTools::sce_to_anndata(sce,
                           anndata_file = opt$output_h5_file)
