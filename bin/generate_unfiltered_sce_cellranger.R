#!/usr/bin/env Rscript

# This script takes the raw H5AD file output by cellranger and
# returns the unfiltered counts matrices as a SingleCellExperiment stored in a .rds file

# import libraries
suppressPackageStartupMessages({
  library(optparse)
  library(SingleCellExperiment)
})
# set up arguments
option_list <- list(
  make_option(
    opt_str = c("--cellranger_file"),
    type = "character",
    help = "path to raw H5AD file output by cellranger"
  ),
  make_option(
    opt_str = c("-u", "--unfiltered_file"),
    type = "character",
    default = "",
    help = "path to output unfiltered rds file. Must end in .rds"
  ),
  make_option(
    opt_str = c("--versions_file"),
    type = "character",
    default = "",
    help = "path to file containing celltype version used for processing"
  ),
  make_option(
    opt_str = c("--reference_index"),
    type = "character",
    default = "",
    help = "path to reference index used with cellranger"
  ),
  make_option(
    opt_str = c("--reference_probeset"),
    type = "character",
    default = "",
    help = "path to reference probe set used withi cellranger"
  ),
  make_option(
    opt_str = c("-t", "--technology"),
    type = "character",
    help = "sequencing technology string to store in metadata"
  ),
  make_option(
    opt_str = c("--assay_ontology_term_id"),
    type = "character",
    default = NULL,
    help = "Experimental Factor Ontology term associated with provided tech_version"
  ),
  make_option(
    opt_str = c("-s", "--seq_unit"),
    type = "character",
    help = "sequencing unit string to store in metadata (e.g., cell, nucleus)"
  ),
  make_option(
    opt_str = c("--library_id"),
    type = "character",
    help = "library id"
  ),
  make_option(
    opt_str = c("--sample_id"),
    type = "character",
    help = "sample id(s). If more than one, separated by commas and/or semicolons."
  ),
  make_option(
    opt_str = c("--project_id"),
    type = "character",
    help = "project id"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

# check that input file exists and output file name ends in rds
stopifnot(
  "--cellranger_file does not exist" = file.exists(opt$cellranger_file),
  "unfiltered file name must end in .rds" = stringr::str_ends(opt$unfiltered_file, ".rds")
)

# parse sample id list
sample_ids <- unlist(stringr::str_split(opt$sample_id, ",|;")) |> sort()

# create sce
unfiltered_sce <- DropletUtils::read10xCounts(
  opt$cellranger_file,
  col.names = TRUE
)

# update the column names to remove the -1 and be consistent with other quantifiers
colnames(unfiltered_sce) <- stringr::str_extract(colnames(unfiltered_sce), "^([ACGT]+)")

# remove existing colData
colData(unfiltered_sce) <- NULL

# get cellranger version
if (file.exists(opt$versions_file)) {
  versions <- jsonlite::fromJSON(opt$versions_file)
  cellranger_version <- versions$pipelines
} else {
  cellranger_version <- NA
}

# make metadata list with scpca information and add to object
metadata_list <- list(
  library_id = opt$library_id,
  sample_id = sample_ids,
  project_id = opt$project_id,
  cellranger_version = cellranger_version,
  reference_index = basename(opt$reference_index),
  reference_probeset = basename(opt$reference_probeset),
  # TODO: Add in total number of reads and mapped reads?
  # Total reads live in metrics_summary.csv and mapped reads is provided as % of reads in cells
  total_reads = NA,
  mapped_reads = NA,
  mapping_tool = "cellranger-multi",
  cellranger_num_cells = nrow(unfiltered_sce),
  tech_version = opt$technology,
  assay_ontology_term_id = opt$assay_ontology_term_id,
  seq_unit = opt$seq_unit,
  transcript_type = "total"
)

metadata(unfiltered_sce) <- metadata_list


# write to rds
readr::write_rds(unfiltered_sce, opt$unfiltered_file, compress = "bz2")
