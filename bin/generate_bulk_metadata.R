#!/usr/bin/env Rscript

# This script generates the metadata.tsv file for bulk RNA-seq libraries.
# It outputs one metadata file per project.

# import libraries
library(optparse)
suppressPackageStartupMessages(library(SingleCellExperiment))

# set up arguments
option_list <- list(
  make_option(
    opt_str = c("-p", "--project_id"),
    type = "character",
    help = "scpca project ID",
  ),
  make_option(
    opt_str = c("-s", "--salmon_dirs"),
    type = "character",
    help = "Path to text file containing salmon output directories, one per line."
  ),
  make_option(
    opt_str = c("--library_metadata_file"),
    type = "character",
    help = "path to metadata file containing scpca_library_id, scpca_sample_id and associated metadata"
  ),
  make_option(
    opt_str = "--metadata_output",
    default = "metadata.tsv",
    type = "character",
    help = "path to metadata tsv output file"
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

# check for project id
if(is.null(opt$project_id)){
  stop("A `project_id` is required.")
}


# replace workflow url and commit if not provided
if (is.null(opt$workflow_url)){
  opt$workflow_url <- NA
}
if (is.null(opt$workflow_version)){
  opt$workflow_version <- NA
}
if (is.null(opt$workflow_commit)){
  opt$workflow_commit <- NA
}

# read in library metadata file 
library_metadata <- readr::read_tsv(opt$library_metadata_file)

# list of paths to salmon log files 
library_ids <- readLines(opt$salmon_dirs)

# subset library metadata to only contain libraries that are being processed 
# only keep metadata columns of interest 
bulk_metadata_df <- library_metadata |>
  dplyr::filter(scpca_library_id %in% library_ids &
                scpca_project_id %in% opt$project_id) |>
  dplyr::select(scpca_sample_id, scpca_library_id, scpca_project_id,
                technology, seq_unit) |> 
  # rename column names to match format of metadata files from other modalities
  dplyr::rename(
    sample_id = scpca_sample_id,
    library_id = scpca_library_id,
    project_id = scpca_project_id
  ) |>
  # add columns with processing information and date processed (same for all libraries )
  dplyr::mutate(date_processed = lubridate::format_ISO8601(lubridate::now(tzone = "UTC"), usetz = TRUE),
                genome_assembly = opt$genome_assembly, 
                workflow = opt$workflow_url,
                workflow_version = opt$workflow_version,
                workflow_commit = opt$workflow_commit)


# add salmon version, index, total_reads, mapped_reads for each library 
add_processing_info <- function(library_id) {
  
  # read in json files for that library containing individual process information
  cmd_info_file <- file.path(library_id, "cmd_info.json")
  meta_info_file <- file.path(library_id, "aux_info", "meta_info.json")
  
  cmd_info <- jsonlite::read_json(cmd_info_file)
  meta_info <- jsonlite::read_json(meta_info_file)
  
  library_processing <- data.frame(
    library_id = library_id,
    salmon_version = cmd_info$salmon_version,
    mapping_index = cmd_info$index,
    total_reads = meta_info$num_processed,
    mapped_reads = meta_info$num_mapped
  )
  
}

bulk_processing_metadata <- purrr::map(library_ids, add_processing_info) |>
  dplyr::bind_rows()

bulk_metadata_df <- bulk_metadata_df |>
  dplyr::left_join(bulk_processing_metadata, by = c("library_id"))

# write out file 
readr::write_tsv(bulk_metadata_df, file = opt$metadata_output)
