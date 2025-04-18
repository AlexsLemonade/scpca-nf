#!/usr/bin/env Rscript

# This script generates a QC report using scpcaTools from a pair of filtered and unfiltered SCE objects

# import libraries
suppressPackageStartupMessages({
  library(optparse)
  library(SingleCellExperiment)
})

# define helper functions #

total_spliced <- function(sce) {
  # get total spliced counts if they exists, otherwise return 0
  if ("spliced" %in% assayNames(sce)) {
    return(sum(assay(sce, "spliced")))
  } else {
    return(0)
  }
}

altexp_totals <- function(sce) {
  # get altExp totals from an SCE as a named list
  exp_totals <- altExpNames(sce) |>
    purrr::set_names() |>
    purrr::map(\(exp) {
      altExp(sce, exp) |>
        counts() |>
        sum()
    })
  return(exp_totals)
}

# set up arguments
option_list <- list(
  make_option(
    opt_str = "--metadata_json",
    type = "character",
    help = "path to metadata json file",
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
    help = "path to output metrics json file",
    default = "metrics.json"
  )
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
  mapped_reads = metadata$mapped_reads,
  unfiltered_cells = metadata$unfiltered_cells,
  filtered_cells = metadata$filtered_cells,
  processed_cells = metadata$processed_cells,
  droplet_filtering_method = metadata$droplet_filtering_method,
  normalization_method = metadata$normalization_method,
  cell_filtering_method = metadata$cell_filtering_method
)

# read sce files and compile metrics for output files
# do these one at a time for memory efficiency
unfiltered_sce <- readr::read_rds(opt$unfiltered_sce)

metrics$unfiltered_total_counts <- sum(counts(unfiltered_sce))
metrics$unfiltered_total_spliced <- total_spliced(unfiltered_sce)
metrics$unfiltered_expressed_genes <- sum(rowSums(counts(unfiltered_sce)) > 0)
metrics$unfiltered_altexp_total <- altexp_totals(unfiltered_sce)

rm(unfiltered_sce)

# add filtered file metrics
if (file.size(opt$filtered_sce) > 0) {
  filtered_sce <- readr::read_rds(opt$filtered_sce)
  metrics$filtered_total_counts <- sum(counts(filtered_sce))
  metrics$filtered_total_spliced <- total_spliced(filtered_sce)
  metrics$filtered_expressed_genes <- sum(rowSums(counts(filtered_sce)) > 0)
  metrics$filtered_altexp_total <- altexp_totals(filtered_sce)
  metrics$miqc_pass_count <- sum(filtered_sce$miQC_pass)
  metrics$scpca_filter_count <- sum(filtered_sce$scpca_filter == "Keep")
  metrics$adt_scpca_filter_count <- sum(filtered_sce$adt_scpca_filter == "Keep")

  rm(filtered_sce)
}

# add processed file metrics
if (file.size(opt$processed_sce) > 0) {
  processed_sce <- readr::read_rds(opt$processed_sce)
  metrics$processed_total_counts <- sum(counts(processed_sce))
  metrics$processed_total_spliced <- total_spliced(processed_sce)
  metrics$processed_total_logcounts <- sum(logcounts(processed_sce))
  metrics$processed_expressed_genes <- sum(rowSums(counts(processed_sce)) > 0)
  metrics$processed_altexp_total <- altexp_totals(processed_sce)
  metrics$hv_genes <- metadata(processed_sce)$highly_variable_genes
  metrics$cluster_algorithm <- metadata(processed_sce)$cluster_algorithm
  # cluster counts as unnamed vector
  metrics$cluster_sizes <- as.vector(table(processed_sce$cluster))
  metrics$singler_reference <- metadata(processed_sce)$singler_reference
  metrics$cellassign_reference <- metadata(processed_sce)$cellassign_reference
  # convert celltype annotation counts to named lists
  metrics$singler_celltypes <- as.list(table(processed_sce$singler_celltype_ontology))
  metrics$cellassign_celltypes <- as.list(table(processed_sce$cellassign_celltype_annotation))
}

jsonlite::write_json(
  metrics,
  path = opt$metrics_json,
  auto_unbox = TRUE,
  pretty = TRUE
)
