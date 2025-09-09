#!/usr/bin/env Rscript

# This script runs adds inferCNV results to an SCE object, including:
# metadata field `infercnv_table`: (data frame) The full table output from the inferCNV HMM
# metadata field `infercnv_options`: (list) The options used when calling inferCNV, obtained
#  from the @options slot of the inferCNV output object
# colData column `total_cnv`: The sum of CNV per cell, calculated from the HMM output

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
    This file should contain a named list with fields `infercnv_table` and `infercnv_options`."
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

# export updated SCE file
readr::write_rds(sce, opts$output_sce_file, compress = "bz2")
