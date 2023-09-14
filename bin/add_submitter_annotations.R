#!/usr/bin/env Rscript

# This script adds submitter annotations, if provided, to the unfiltered SCE's colData.
# If the file was not provided, this script has no effect.

# import libraries
suppressPackageStartupMessages({
  library(optparse)
  library(SingleCellExperiment)
})
# set up arguments
option_list <- list(
  make_option(
    opt_str = c("-f", "--sce_file"),
    type = "character",
    help = "path to SingleCellExperiment file to update. Must end in .rds"
  ),
  make_option(
    opt_str = c("--library_id"),
    type = "character",
    help = "library id"
  ),
  make_option(
    opt_str = c("--submitter_cell_types_file"),
    type = "character",
    help = "path to tsv file containing submitter-provided cell type annotations"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

# check that output file name ends in .rds
if (!(stringr::str_ends(opt$unfiltered_file, ".rds"))) {
  stop("unfiltered file name must end in .rds")
}


# check that submitter cell types file exists
if (!file.exists(opt$submitter_cell_types_file)) {
  stop("submitter cell type annotations file not found.")
}


# Read in unfiltered sce
sce <- readr::read_rds(opt$unfiltered_file)

# Read in celltypes, and filter to relevant information
submitter_cell_types_df <- readr::read_tsv(opt$submitter_cell_types_file) |>
  dplyr::filter(scpca_library_id == opt$library_id) |>
  dplyr::select(
    barcodes = cell_barcode,
    submitter_celltype_annotations = cell_type_assignment
  ) |>
   # in the event of NA values, change to "Unclassified cell"
    tidyr::replace_na(list(submitter_celltype_annotations = "Unclassified cell"))
  )
  
# Join in submitter cell types
colData(sce) <- colData(sce) |>
  as.data.frame() |>
  dplyr::left_join(
    submitter_cell_types_df
  ) |>
  # make sure we keep rownames
  DataFrame(row.names = colData(sce)$barcodes)

# Indicate that we have submitter celltypes in metadata, 
#  saving in same spot as for actual celltyping
metadata(sce)$celltype_methods <- c(metadata(sce)$celltype_methods, "submitter")

# Write unfiltered RDS back to file
readr::write_rds(sce, opt$unfiltered_file, compress = "gz")
