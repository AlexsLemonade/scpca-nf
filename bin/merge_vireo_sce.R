#!/usr/bin/env Rscript

# import libraries
suppressPackageStartupMessages({
  library(optparse)
  library(SingleCellExperiment)
})
# set up arguments
option_list <- list(
  make_option(
    opt_str = c("-s", "--sce_file"),
    type = "character",
    help = "path to rds file with an sce object"
  ),
  make_option(
    opt_str = c("-o", "--output_sce_file"),
    type = "character",
    help = "path to output with multiplex data added. Must end in .rds"
  ),
  make_option(
    opt_str = c("-v", "--vireo_dir"),
    type = "character",
    help = "path to vireo output directory",
    default = NULL
  )
)

opt <- parse_args(OptionParser(option_list = option_list))


# check that unfiltered file file exists
if(!file.exists(opt$sce_file)){
  stop("Missing unfiltered.rds file")
}

# check that output file name ends in .rds
if(!(stringr::str_ends(opt$output_sce_file, ".rds"))){
  stop("filtered file name must end in .rds")
}

# check for donr_ids file in vireo_dir
vireo_file <- file.path(opt$vireo_dir, "donor_ids.tsv")
if(!file.exists(vireo_file)){
  stop("Missing donor_ids.tsv file in vireo directory")
}

# read in sce rds file
sce <- readr::read_rds(opt$sce_file)

# read in vireo file
vireo_table <- readr::read_tsv(vireo_file)

# add vireo results to sce
sce <- scpcaTools::add_demux_vireo(sce, vireo_table)

# write filtered sce to output
readr::write_rds(sce, opt$output_sce_file, compress = "gz")
