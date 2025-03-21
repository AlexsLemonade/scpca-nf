#!/usr/bin/env Rscript

# This script generates a QC report using scpcaTools from a pair of filtered and unfiltered SCE objects

# import libraries
suppressPackageStartupMessages({
  library(optparse)
  library(SingleCellExperiment)
})

# set up arguments
option_list <- list(
  make_option(
    opt_str = "--metadata_json",
    default = "metadata.json",
    type = "character",
    help = "path to metadata json  file",
    default = ""
  ),
  make_option(
    opt_str = c("-u", "--unfiltered_sce"),
    type = "character",
    help = "path to rds file with unfiltered sce object",
    default = ""
  ),
  make_option(
    opt_str = c("-f", "--filtered_sce"),
    type = "character",
    help = "path to rds file with filtered sce object",
    default = ""
  ),
  make_option(
    opt_str = c("-p", "--processed_sce"),
    type = "character",
    help = "path to rds file with processed sce object",
    default = ""
  ),
  make_option(
    opt_str = c("-o", "--metrics_json"),
    type = "character",
    help = "path to output metrics json file"
  ),
)

opt <- parse_args(OptionParser(option_list = option_list))

stopifnot(
  "Unfiltered .rds file missing or `unfiltered_sce` not specified." = file.exists(opt$unfiltered_sce),
  "metadata json file missing or `metadata_json` not specified." = file.exists(opt$metadata_json)
)


# Read in metadata.json
metadata <- jsonlite::read_json(opt$metadata_json)

# create initial metrics object
# populate with metrics from metadata
metrics <- list(
  project_id = metadata$project_id,
  library_id = metadata$library_id,
  sample_id = metadata$sample_id,
  unfiltered_cells = metadata$unfiltered_cells,
  processed_cells = metadata$processed_cells,
  mapped_reads = metadata$mapped_reads,
  droplet_filtering_method = metadata$droplet_filtering_method,
  normalization_method = metadata$normalization_method,
  cell_filtering_method = metadata$cell_filtering_method
)

# read sce files and compile metrics for output files
# do these one at a time for memory efficiency
unfiltered_sce <- readr::read_rds(opt$unfiltered_sce)

metrics$unfiltered_total_counts <- counts(unfiltered_sce) |> sum()
metrics$unfiltered_expressed_genes <- rowsum(counts(unfiltered_sce)) > 0 |> sum()

rm(unfiltered_sce)

# make sure filtered sce has an object, otherwise set to NULL
if (file.size(opt$filtered_sce) > 0) {
  filtered_sce <- readr::read_rds(opt$filtered_sce)
  metrics$filtered_total_counts <- counts(filtered_sce) |> sum()
  metrics$filtered_expressed_genes <- rowsum(counts(filtered_sce)) > 0 |> sum()
  metrics$miqc_pass_count <- sum(filtered_sce$miqc_pass)

  rm(filtered_sce)
}

# make sure processed sce has an object, otherwise set to NULL
if (file.size(opt$processed_sce) > 0) {
  processed_sce <- readr::read_rds(opt$processed_sce)
  processed_sce_meta <- metadata(processed_sce)
} else {
  processed_sce <- NULL
  processed_sce_meta <- NULL
}
