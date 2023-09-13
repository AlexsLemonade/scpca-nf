#!/usr/bin/env Rscript

# This script generates a QC report using scpcaTools from a pair of filtered and unfiltered SCE objects

# import libraries
library(optparse)
suppressPackageStartupMessages(library(SingleCellExperiment))

# set up arguments
option_list <- list(
  make_option(
    opt_str = c("--report_template"),
    type = "character",
    default = NULL,
    help = "path to rmd template file for report"
  ),
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
    opt_str = c("-p", "--processed_sce"),
    type = "character",
    help = "path to rds file with processed sce object"
  ),
  make_option(
    opt_str = c("-l", "--library_id"),
    type = "character",
    help = "Library identifier for report and metadata file"
  ),
  make_option(
    opt_str = c("-s", "--sample_id"),
    type = "character",
    help = "Sample identifier(s) for metadata file. For a multiplexed library, a comma or semicolon separated list."
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
  ),
  make_option(
    opt_str = "--demux_method",
    type = "character",
    default = "vireo",
    help = "Demultiplexing method to use for multiplexed samples. One of `vireo`, `HTOdemux`, or `HashedDrops`"
  ),
  make_option(
    opt_str = "--seed",
    type = "integer",
    default = NULL,
    help = "Optional random seed used in certain QC report visualizations."
  )
)

opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$library_id)) {
  stop("A `library_id` is required.")
}

# check that template file, if given, exists
if (!is.null(opt$report_template) && !file.exists(opt$report_template)) {
  stop("Specified `report_template` could not be found.")
}

if (is.null(opt$unfiltered_sce) || !file.exists(opt$unfiltered_sce)) {
  stop("Unfiltered .rds file missing or `unfiltered_sce` not specified.")
}
if (is.null(opt$filtered_sce) || !file.exists(opt$filtered_sce)) {
  stop("Filtered .rds file missing or `filtered_sce` not specified.")
}
demux_methods <- c("vireo", "HTODemux", "HashedDrops")
if (!opt$demux_method %in% demux_methods) {
  stop("Unknown `demux_method` value. Must be one of `vireo`, `HTOdemux`, or `HashedDrops`")
}

if (opt$workflow_url == "null") {
  opt$workflow_url <- NA
}
if (opt$workflow_version == "null") {
  opt$workflow_version <- NA
}
if (opt$workflow_commit == "null") {
  opt$workflow_commit <- NA
}

# read sce files
unfiltered_sce <- readr::read_rds(opt$unfiltered_sce)
filtered_sce <- readr::read_rds(opt$filtered_sce)
processed_sce <- readr::read_rds(opt$processed_sce)

# Compile metadata for output files
sce_meta <- metadata(unfiltered_sce)
filtered_sce_meta <- metadata(filtered_sce)
processed_sce_meta <- metadata(processed_sce)

# Parse sample ids
sample_ids <- unlist(stringr::str_split(opt$sample_id, ",|;")) |> sort()

# check for multiplexing
multiplexed <- if (length(sample_ids) > 1) {
  TRUE
} else {
  FALSE
}

# sanity check ids
if (!is.null(sce_meta$sample_id)) {
  if (!all.equal(sample_ids, sce_meta$sample_id)) {
    stop("--sample_id  does not match SCE metadata")
  }
}
if (!is.null(sce_meta$library_id)) {
  if (opt$library_id != sce_meta$library_id) {
    stop("--library_id  does not match SCE metadata")
  }
}

# check for alt experiments (CITE-seq, etc)
alt_expts <- altExpNames(unfiltered_sce)
has_citeseq <- "adt" %in% alt_expts
has_cellhash <- "cellhash" %in% alt_expts


metadata_list <- list(
  library_id = opt$library_id,
  sample_id = opt$sample_id,
  technology = opt$technology,
  seq_unit = opt$seq_unit,
  is_multiplexed = multiplexed,
  has_citeseq = has_citeseq,
  has_cellhash = has_cellhash,
  filtered_cells = ncol(filtered_sce),
  unfiltered_cells = ncol(unfiltered_sce),
  droplet_filtering_method = filtered_sce_meta$filtering_method,
  total_reads = sce_meta$total_reads,
  mapped_reads = sce_meta$mapped_reads,
  genome_assembly = opt$genome_assembly,
  mapping_index = sce_meta$reference_index,
  transcript_type = sce_meta$transcript_type,
  cell_filtering_method = processed_sce_meta$scpca_filter_method,
  normalization_method = processed_sce_meta$normalization,
  min_gene_cutoff = processed_sce_meta$min_gene_cutoff,
  prob_compromised_cutoff = processed_sce_meta$prob_compromised_cutoff,
  date_processed = lubridate::format_ISO8601(lubridate::now(tzone = "UTC"), usetz = TRUE),
  salmon_version = sce_meta$salmon_version,
  alevin_fry_version = sce_meta$alevinfry_version,
  workflow = opt$workflow_url,
  workflow_version = opt$workflow_version,
  workflow_commit = opt$workflow_commit
) |>
  purrr::map(\(x) {
    if (is.null(x)) NA else x
  }) # convert any NULLS to NA

# add adt methods if citeseq
if (has_citeseq) {
  metadata_list <- append(
    metadata_list,
    list(
      adt_filtering_method = processed_sce_meta$adt_scpca_filter_method,
      adt_normalization_method = processed_sce_meta$adt_normalization
    )
  )
}


# estimate cell counts for multiplexed samples
if (multiplexed) {
  demux_column <- paste0(opt$demux_method, "_sampleid")
  demux_counts <- colData(filtered_sce)[[demux_column]] |>
    table() |>
    as.list() # save as a list for json output

  # add demux info to the metadata list
  metadata_list <- append(
    metadata_list,
    list(
      demux_method = opt$demux_method,
      demux_samples = sample_ids,
      sample_cell_estimates = demux_counts
    )
  )
}

# Output metadata as JSON
jsonlite::write_json(metadata_list, path = opt$metadata_json, auto_unbox = TRUE, pretty = TRUE)

scpcaTools::generate_qc_report(
  library_id = metadata_list$library_id,
  unfiltered_sce = unfiltered_sce,
  filtered_sce = filtered_sce,
  processed_sce = processed_sce,
  report_template = opt$report_template,
  output = opt$qc_report_file,
  extra_params = list(seed = opt$seed)
)
