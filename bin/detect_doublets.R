#!/usr/bin/env Rscript

# Script used to perform doublet detection using scDblFinder on a given SCE object
#
# This script reads in an RDS file containing an SCE and runs scDblFinder.
# The full results are stored in the SCE's metadata, and classifications are
# stored in the SCE's colData slot.
# The SCE is then exported.

# import libraries
library(optparse)
library(SingleCellExperiment)

# Set up optparse options
option_list <- list(
  make_option(
    opt_str = c("--input_sce_file"),
    type = "character",
    default = "",
    help = "Path to RDS file that contains the SCE object to perform doublet detection on."
  ),
  make_option(
    opt_str = c("-o", "--output_sce_file"),
    type = "character",
    default = "",
    help = "Output path for updated SCE file. Must end in .rds"
  ),
  make_option(
    opt_str = c("-t", "--threads"),
    type = "integer",
    default = 4,
    help = "Number of multiprocessing threads to use."
  ),
  make_option(
    opt_str = c("--random_seed"),
    type = "integer",
    default = 2024,
    help = "Random seed for reproducibility"
  )
)

# Setup ------------------------------
opt <- parse_args(OptionParser(option_list = option_list))
set.seed(opt$random_seed)

# check and read in SCE file
if (!file.exists(opt$input_sce_file)) {
  stop("Input `input_sce_file` is missing.")
}
sce <- readr::read_rds(opt$input_sce_file)

if (!stringr::str_ends(opt$output_sce_file, ".rds")) {
  stop("`output_sce_file` must end in .rds")
}

# set up multiprocessing params
if (opt$threads > 1) {
  bp_param <- BiocParallel::MulticoreParam(opt$threads)
} else {
  bp_param <- BiocParallel::SerialParam()
}

# run scDblFinder -----------
doublet_result <- scDblFinder::scDblFinder(
  sce,
  BPPARAM = bp_param,
  returnType = "table"
)

# store the full table result in metadata
metadata(sce)$scDblFinder_result <- doublet_result

# store the `score` and `class` columns in the colData
doublet_result_cells <- doublet_result |>
  as.data.frame() |>
  tibble::rownames_to_column("barcodes") |>
  # keep only actual barcodes to remove the artificial doublets
  dplyr::filter(barcodes %in% colnames(sce)) |>
  dplyr::select(barcodes, scDblFinder_class = class, scDblFinder_score = score)

colData(sce) <- colData(sce) |>
  as.data.frame() |>
  dplyr::left_join(doublet_result_cells, by = "barcodes") |>
  DataFrame(row.names = colnames(sce))

# Export updated SCE with doublet information
readr::write_rds(sce, opt$output_sce_file, compress = "bz2")
