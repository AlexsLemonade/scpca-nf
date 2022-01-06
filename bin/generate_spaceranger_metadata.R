#!/usr/bin/env Rscript

# This script generates the metadata.json file for Spatial Transcriptomics libraries. 

# import libraries
library(optparse)
suppressPackageStartupMessages(library(SingleCellExperiment))

# set up arguments
option_list <- list(
  make_option(
    opt_str = c("-l", "--library_id"),
    type = "character",
    help = "Library identifier for report and metadata file"
  ),
  make_option(
    opt_str = c("-s", "--sample_id"),
    type = "character",
    help = "Sample identifier for metadata file"
  ),
  make_option(
    opt_str = c("-u", "--unfiltered_barcodes_file"),
    type = "character",
    help = "path to file containing list of all spot barcodes identified in library with no filtering"
  ),
  make_option(
    opt_str = c("-f", "--filtered_barcodes_file"),
    type = "character",
    help = "path to file containing list of spot barcodes identified in library after filtering"
  ),
  make_option(
    opt_str = c("--metrics_summary_file"),
    type = "character",
    help = "path to file containing summary of run metrics from spaceranger count, the metrics_summary.csv file output by spaceranger count"
  ),
  make_option(
    opt_str = c("--spaceranger_versions_file"),
    type = "character",
    help = "path to json file containing spaceranger version information, the _versions file output by spaceranger count"
  ),
  make_option(
    opt_str = "--metadata_json",
    default = "metadata.json",
    type = "character",
    help = "path to metadata json output file"
  ),
  make_option(
    opt_str = "--technology",
    type = "character",
    default = NA,
    help = "sequencing technology"
  ),
  make_option(
    opt_str = "--seq_unit",
    type = "character",
    default = NA,
    help = "sequencing unit"
  ),
  make_option(
    opt_str = "--genome_assembly",
    type = "character",
    default = NA,
    help = "genome assembly used for mapping"
  ),
  make_option(
    opt_str = "--index_path",
    type = "character",
    default = NA,
    help = "path to index used for mapping"
  ),
  make_option(
    opt_str = "--workflow_url",
    type = "character",
    default = NA,
    help = "workflow github url"
  ),
  make_option(
    opt_str = "--workflow_version",
    type = "character",
    default = NA,
    help = "workflow version identifier"
  ),
  make_option(
    opt_str = "--workflow_commit",
    type = "character",
    default = NA,
    help = "workflow commit hash"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))
if(is.null(opt$library_id)){
  stop("A `library_id` is required.")
}
# check that barcode files exist
if(is.null(opt$unfiltered_barcodes_file) || !dir.exists(opt$unfiltered_barcodes_file)){
  stop("Unfiltered barcodes file missing or `unfiltered_barcodes_file` not specified.")
}
if(is.null(opt$filtered_barcodes_file) || !file.exists(opt$filtered_barcodes_file)){
  stop("Unfiltered barcodes file missing or `unfiltered_barcodes_file` not specified.")
}

# check that metrics summary file exists 
if(is.null(opt$metrics_summary_file) || !file.exists(opt$metrics_summary_file)){
  stop("Metrics summary file missing or `metrics_summary_file` not specified.")
}

# check that version file exists
if(is.null(opt$spaceranger_versions_file) || !file.exists(opt$spaceranger_versions_file)){
  stop("Versions file missing or `spaceranger_versions_file` not specified.")
}

if (opt$workflow_url == "null"){
  opt$workflow_url <- NA
}
if (opt$workflow_version == "null"){
  opt$workflow_version <- NA
}
if (opt$workflow_commit == "null"){
  opt$workflow_commit <- NA
}

# read in barcode files
unfiltered_barcodes <- readr::read_tsv(opt$unfiltered_barcodes_file,
                                       col_names = c("barcode"))
filtered_barcodes <- readr::read_tsv(opt$filtered_barcodes_file,
                                     col_names = c("barcode"))

# read in metrics summary 
metrics_summary <- readr::read_csv(opt$metrics_summary_file)

# read in versions file 
spaceranger_versions <- jsonlite::read_json(opt$spaceranger_versions_file)

# obtain just the index name and drop the rest of the path 
index_name <- basename(opt$index_path)

# compile metadata list 
metadata_list <- list(
  library_id = opt$library_id,
  sample_id = opt$sample_id,
  technology = opt$technology,
  seq_unit = opt$seq_unit,
  filtered_spots = nrow(filtered_barcodes),
  unfiltered_spots = nrow(unfiltered_barcodes),
  total_reads = metrics_summary$`Number of Reads`,
  mapped_reads = metrics_summary$`Reads Mapped Confidently to Genome`,
  tissue_spots = metrics_summary$`Number of Spots Under Tissue`,
  genome_assembly = opt$genome_assembly,
  mapping_index = index_name,
  date_processed = lubridate::format_ISO8601(lubridate::now(tzone = "UTC"), usetz = TRUE),
  spaceranger_version = spaceranger_versions$pipelines,
  workflow = opt$workflow_url,
  workflow_version = opt$workflow_version,
  workflow_commit = opt$workflow_commit
) |>
  purrr::map(~if(is.null(.)) NA else .) # convert any NULLS to NA

# Output metadata as JSON
jsonlite::write_json(metadata_list, path = opt$metadata_json, auto_unbox = TRUE)

