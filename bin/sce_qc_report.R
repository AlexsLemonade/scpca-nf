#!/usr/bin/env Rscript

# This script generates a QC report using scpcaTools from a pair of filtered and unfiltered SCE objects

# import libraries
library(optparse)
suppressPackageStartupMessages(library(SingleCellExperiment))

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
    opt_str = c("-q", "--qc_report_file"),
    default = "qc_report.html",
    type = "character",
    help = "path to QC report output file"
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
# check that input files exist
if(is.null(opt$unfiltered_sce) || !file.exists(opt$unfiltered_sce)){
  stop("Unfiltered .rds file missing or `unfiltered_sce` not specified.")
}
if(is.null(opt$filtered_sce) || !file.exists(opt$filtered_sce)){
  stop("Filtered .rds file missing or `filtered_sce` not specified.")
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

# read sce files
unfiltered_sce <- readr::read_rds(opt$unfiltered_sce)
filtered_sce <- readr::read_rds(opt$filtered_sce)

# Compile metadata for output files
sce_meta <- metadata(unfiltered_sce)
filtered_sce_meta <- metadata(filtered_sce)

# check for alt experiments (CITE-seq, etc)
alt_expts <- altExpNames(unfiltered_sce)
has_citeseq <- "CITEseq" %in% alt_expts

metadata_list <- list(
  library_id = opt$library_id,
  sample_id = opt$sample_id,
  technology = opt$technology,
  seq_unit = opt$seq_unit,
  has_citeseq = has_citeseq,
  filtered_cells = ncol(filtered_sce),
  unfiltered_cells = ncol(unfiltered_sce),
  total_reads = sce_meta$total_reads,
  mapped_reads = sce_meta$mapped_reads,
  filtering_method = filtered_sce_meta$filtering_method,
  genome_assembly = opt$genome_assembly,
  mapping_index = sce_meta$reference_index,
  transcript_type = sce_meta$transcript_type,
  date_processed = lubridate::format_ISO8601(lubridate::now(tzone = "UTC"), usetz = TRUE),
  salmon_version = sce_meta$salmon_version,
  alevin_fry_version = sce_meta$alevinfry_version,
  workflow = opt$workflow_url,
  workflow_version = opt$workflow_version,
  workflow_commit = opt$workflow_commit
) |>
  purrr::map(~if(is.null(.)) NA else .) # convert any NULLS to NA

# Output metadata as JSON
jsonlite::write_json(metadata_list, path = opt$metadata_json, auto_unbox = TRUE)
  
scpcaTools::generate_qc_report(
  sample_name = metadata_list$library_id,
  unfiltered_sce = unfiltered_sce,
  filtered_sce = filtered_sce,
  output = opt$qc_report_file
)

