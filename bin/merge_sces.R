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

# list of paths to salmon files 
input_sce_files <- stringr::str_split(opt$infiles, ',')
# pull library ids from list of 
input_library_ids <- stringr::str_split(opt$library_ids, ',')

# get list of sces
sce_list <- purrr::map(input_sce_files, readr::read_rds)
names(sce_list) <- input_library_ids

merged_sce <- sce_list[[1]]

# write out merged sce file 
readr::write_rds(merged_sce, opt$output_sce_file)
