#!/usr/bin/env Rscript

# This script adds inferCNV results to an SCE object, including:
# metadata field `infercnv_table`: (data frame) The full table output from the inferCNV HMM
# metadata field `infercnv_options`: (list) The options used when calling inferCNV, obtained
#  from the @options slot of the inferCNV output object
# colData column `total_cnv`: The sum of CNV per cell, calculated from the HMM output
#
# If the inputted inferCNV results file is empty, instead only an info message
# for why there are no results is added to the SCE metadata field `infercnv_options`

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(optparse)
})

option_list <- list(
  make_option(
    opt_str = "--input_sce_file",
    type = "character",
    default = "",
    help = "Path to the SCE file to run inferCNV on"
  ),
  make_option(
    opt_str = "--infercnv_results_file",
    type = "character",
    default = "",
    help = "Path to the inferCNV results file.
    If not empty, this file should contain a named list with fields `infercnv_table` and `infercnv_options`."
  ),
  make_option(
    opt_str = "--infercnv_threshold",
    type = "integer",
    help = "Minimum number of cells needed to run inferCNV, used here to report null results."
  ),
  make_option(
    opt_str = "--output_sce_file",
    type = "character",
    help = "Path to the output SCE with inferCNV information included"
  )
)

opts <- parse_args(OptionParser(option_list = option_list))

stopifnot(
  "input_sce_file does not exist" = file.exists(opts$input_sce_file),
  "infercnv_results_file does not exist" = file.exists(opts$infercnv_results_file),
  "output_sce_file was not provided" = !is.null(opts$output_sce_file)
)

# read SCE file
sce <- readRDS(opts$input_sce_file)

# check if we have inferCNV results based on file size
if (file.info(opts$infercnv_results_file)$size > 0) {
  # read inferCNV results
  infercnv_results <- readRDS(opts$infercnv_results_file)

  # add information to metadata
  metadata(sce)$infercnv_options <- infercnv_results$infercnv_options
  metadata(sce)$infercnv_table <- infercnv_results$infercnv_table

  # add a total_cnv column to colData
  total_cnv_df <- infercnv_results$infercnv_table |>
    tidyr::pivot_longer(
      starts_with("has_cnv_"),
      names_to = "chr",
      values_to = "cnv"
    ) |>
    dplyr::group_by(barcodes) |>
    dplyr::summarize(total_cnv = sum(cnv))
  colData(sce) <- colData(sce) |>
    as.data.frame() |>
    dplyr::left_join(total_cnv_df, by = "barcodes") |>
    DataFrame(row.names = colnames(sce))
} else {
  # no results; add info message to metadata
  if (sum(sce$is_infercnv_reference) < opts$infercnv_threshold) {
    metadata(sce)$infercnv_options <- glue::glue(
      "inferCNV not run; insufficient cells for a normal reference"
    )
  } else {
    metadata(sce)$infercnv_options <- "inferCNV failed to run"
  }
}

# export updated SCE file
readr::write_rds(sce, opts$output_sce_file, compress = "bz2")
