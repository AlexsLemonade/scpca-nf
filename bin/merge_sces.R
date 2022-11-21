#!/usr/bin/env Rscript

# import libraries
library(magrittr)
library(optparse)
library(SingleCellExperiment)

# Set up optparse options
option_list <- list(
  make_option(
    opt_str = c("--input_library_ids"),
    type = "character",
    help = "Comma separated list of library IDs corresponding to the libraries being integrated."
  ),
  make_option(
    opt_str = c("--input_sce_files"),
    type = "character",
    help = "Comma separated list of input sce file paths corresponding to the sces being integrated."
  ),
  make_option(
    opt_str = c("-o", "--output_sce_file"),
    type = "character",
    help = "Path to output RDS file, must end in .rds"
  )
)

# Setup ------------------------------------------------------------------------

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

# check that file extension for output file is correct 

# list of paths to salmon files 
input_sce_files <- stringr::str_split(opt$input_sce_files, ',')
# pull library ids from list of 
input_library_ids <- stringr::str_split(opt$library_ids, ',')

# check that input files exist 

# get list of sces
sce_list <- purrr::map(input_sce_files, readr::read_rds)
names(sce_list) <- input_library_ids

# create combined SCE object with function in scpcaTools 

# HVG calculation (from perform_hvg_selection in integration-helpers)

# Dim Reduction PCA and UMAP 

# write out merged sce file 
readr::write_rds(sce_list, opt$output_sce_file)
