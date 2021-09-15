#!/usr/bin/env Rscript

# This script generates a QC report using scpcaTools from a pair of filtered and unfiltered SCE objects

# import libraries
library(magrittr)
library(optparse)

# set up arguments
option_list <- list(
  make_option(
    opt_str = c("-u", "--unfiltered_sce"),
    type = "character",
    help = "path to rds file with unfiltered sce object"
  ),
  make_option(
    opt_str = c("-f", "--filtered_sce"),
    type = "character",
    help = "path to rds file with filtered sce object"
  ),
  make_option(
    opt_str = c("-s", "--sample_id"),
    type = "character",
    help = "Sample identifier for report"
  ),
  make_option(
    opt_str = c("-o", "--output_file"),
    default = "qc_report.html",
    type = "character",
    help = "path to output file"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

# check that input files exists
if(!file.exists(opt$unfiltered_sce)){
  stop("Missing unfiltered .rds file")
}
if(!file.exists(opt$filtered_sce)){
  stop("Missing filtered .rds file")
}

# read sce files
unfiltered_sce <- readr::read_rds(opt$unfiltered_sce)
filtered_sce <- readr::read_rds(opt$filtered_sce)

scpcaTools::generate_qc_report(
  sample_name = opt$sample_name,
  unfiltered_sce = unfiltered_sce,
  filtered_sce = filtered_sce,
  output = opt$output_file
)
