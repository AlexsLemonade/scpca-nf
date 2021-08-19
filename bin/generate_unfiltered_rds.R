#!/usr/bin/env Rscript

# This script takes the output folder from alevin-fry as input and
# returns the unfiltered counts matrices as a SingleCellExperiment stored in a .rds file

# import libraries
library(magrittr)
library(optparse)
library(scpcaTools)

# set up arguments
option_list <- list(
  make_option(
    opt_str = c("-s", "--seq_unit"),
    type = "character",
    help = "include counts for spliced cDNA only (spliced) or unspliced and spliced cDNA (unspliced)"
  ),
  make_option(
    opt_str = c("-a", "--alevin_dir"),
    type = "character",
    help = "directory with alevin output files"
  ),
  make_option(
    opt_str = c("-u", "--unfiltered_file"),
    type = "character",
    help = "path to output unfiltered rds file. Must end in .rds"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

# check for compatible sequencing unit types
if(!(opt$seq_unit %in% c("cell", "nucleus"))){
  stop("Sequencing unit must be of type cell or nucleus")
}

# check that output file name ends in .rds
if(!(stringr::str_ends(opt$unfiltered_file, ".rds"))){
  stop("unfiltered file name must end in .rds")
}

# convert seq_unit to spliced or unspliced to determine which types of transcripts to include in final counts matrix
which_counts <- dplyr::case_when(opt$seq_unit == "cell" ~ "spliced",
                                 opt$seq_unit == "nucleus" ~ "unspliced")

# get unfiltered sce
unfiltered_sce <- import_quant_data(quant_dir = opt$alevin_dir,
                                    tool = "alevin-fry",
                                    which_counts = which_counts,
                                    usa_mode = TRUE)

# write to rds
readr::write_rds(unfiltered_sce, opt$unfiltered_file)

