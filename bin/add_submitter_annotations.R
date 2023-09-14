#!/usr/bin/env Rscript

# This script adds submitter annotations, if provided, to the an SCE object's colData.

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
if (!(stringr::str_ends(opt$sce_file, ".rds"))) {
  stop("SingleCellExperiment file name must end in .rds")
}


# check that submitter cell types file exists
if (!file.exists(opt$submitter_cell_types_file)) {
  stop("submitter cell type annotations file not found.")
}


# Read in sce
sce <- readr::read_rds(opt$sce_file)

# Read in celltypes
submitter_cell_types_df <- readr::read_tsv(opt$submitter_cell_types_file) |>
  # filter to relevant information
  dplyr::filter(scpca_library_id == opt$library_id) |>
  dplyr::select(
    barcodes = cell_barcode,
    submitter_celltype_annotations = cell_type_assignment
  ) |>
  # in the event of NA values, change to "Unclassified cell"
  tidyr::replace_na(
    list(submitter_celltype_annotations = "Unclassified cell")
  ) |>
  # join with colData
  dplyr::right_join(
    colData(sce) |>
      as.data.frame()
  )
  
# Check rows before sending back into the SCE object
if (nrow(submitter_cell_types_df) != ncol(sce)) {
  stop("Could not add submitter annotations to SCE object. There should only be one annotation per cell.")
}

# Rejoin with colData, making sure we keep rownames
colData(sce) <- DataFrame(
  submitter_cell_types_df, 
  row.names = colData(sce)$barcodes
)

# Indicate that we have submitter celltypes in metadata, 
#  saving in same spot as for actual celltyping
metadata(sce)$celltype_methods <- c(metadata(sce)$celltype_methods, "submitter")

# Write SCE back to file
readr::write_rds(sce, opt$sce_file, compress = "gz")
