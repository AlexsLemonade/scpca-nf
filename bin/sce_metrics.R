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
  cell_filtering_method = metadata$cell_filtering_method,
  workflow_version = metadata$workflow_version,
  date_processed = metadata$date_processed,
  salmon_version = metadata$salmon_version,
  alevinfry_version = metadata$alevinfry_version,
  min_gene_cutoff = metadata$min_gene_cutoff,
  prob_compromised_cutoff = metadata$prob_compromised_cutoff
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
  metrics$miqc_pass_count <- ifelse(
    is.null(filtered_sce$miQC_pass),
    NA_integer_,
    sum(filtered_sce$miQC_pass)
  )
  metrics$scpca_filter_count <- sum(filtered_sce$scpca_filter == "Keep")
  metrics$adt_scpca_filter_count <- ifelse(
    is.null(filtered_sce$adt_scpca_filter),
    NA_integer_,
    sum(filtered_sce$adt_scpca_filter == "Keep")
  )
  metrics$scdblfinder_total_doublets <- sum(filtered_sce$scDblFinder_class == "doublet")

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
  metrics$infercnv_total_cnv <- sum(processed_sce$infercnv_total_cnv)
  # cluster counts as unnamed vector
  metrics$cluster_sizes <- as.vector(table(processed_sce$cluster))

  # infercnv reference cells as named list
  metrics$infercnv_reference_cells <- colData(processed_sce) |>
    as.data.frame() |>
    dplyr::filter(is_infercnv_reference) |>
    dplyr::pull(consensus_celltype_annotation) |>
    table() |>
    as.list()

  metrics$singler_reference <- ifelse(
    is.null(metadata(processed_sce)$singler_reference),
    NA_character_,
    metadata(processed_sce)$singler_reference
  )
  metrics$cellassign_reference <- ifelse(
    is.null(metadata(processed_sce)$cellassign_reference),
    NA_character_,
    metadata(processed_sce)$cellassign_reference
  )
  metrics$scimilarity_model <- ifelse(
    is.null(metadata(processed_sce)$scimilarity_model),
    NA_character_,
    metadata(processed_sce)$scimilarity_model
  )
  # convert celltype annotation counts to named lists
  metrics$singler_celltypes <- as.list(table(processed_sce$singler_celltype_annotation))
  metrics$cellassign_celltypes <- as.list(table(processed_sce$cellassign_celltype_annotation))
  metrics$scimilarity_celltypes <- as.list(table(processed_sce$scimilarity_celltype_annotation))
  metrics$consensus_celltypes <- as.list(table(processed_sce$consensus_celltype_annotation))
}

jsonlite::write_json(
  metrics,
  path = opt$metrics_json,
  pretty = TRUE
)
