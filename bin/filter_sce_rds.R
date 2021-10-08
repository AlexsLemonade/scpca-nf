#!/usr/bin/env Rscript

# This script takes a SingleCellExperiment stored in a .rds file and
# filters it using emptyDrops, adding miQC metrics for probability compromised

# import libraries
suppressPackageStartupMessages({
  library(optparse)
  library(SingleCellExperiment)
})
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
    opt_str = c("-l", "--lower"),
    type = "integer",
    help = "Value specifying the lower bound on total UMI count used in filtering with DropletUtils::emptyDrops."
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
filtered_sce <- scpcaTools::filter_counts(unfiltered_sce,
                                          lower = opt$lower)
# remove unfiltered for memory saving
rm(unfiltered_sce)

# need to remove old gene-level rowData statistics first
drop_cols = colnames(rowData(filtered_sce)) %in% c('mean', 'detected')
rowData(filtered_sce) <- rowData(filtered_sce)[!drop_cols] 

# recalculate rowData and add to filtered sce
filtered_sce <- filtered_sce |>
  scuttle::addPerFeatureQCMetrics()

# add prob_compromised and miQC model to metadata
filtered_sce <- add_miQC()

# grab names of altExp, if any
alt_names <- altExpNames(filtered_sce)

for (alt in alt_names) {
  # remove old row data from unfiltered
  drop_cols = colnames(rowData(altExp(filtered_sce, alt))) %in% c('mean', 'detected')
  rowData(altExp(filtered_sce, alt)) <- rowData(altExp(filtered_sce, alt))[!drop_cols] 

  # add alt experiment features stats for filtered data
  altExp(filtered_sce, alt) <- scuttle::addPerFeatureQCMetrics(altExp(filtered_sce, alt))
}

# write filtered sce to output
readr::write_rds(filtered_sce, opt$filtered_file, compress = "gz")

