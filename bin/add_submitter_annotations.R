#!/usr/bin/env Rscript

# This script adds submitter annotations, if provided, to an SCE object's colData.

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


# check that submitter cell types file exists and is a TSV
if (!file.exists(opt$submitter_cell_types_file)) {
  stop("Submitter cell type annotations file not found.")
}
if (!(stringr::str_ends(opt$submitter_cell_types_file, ".tsv"))) {
  stop("Submitter cell type annotations file must be a TSV.")
}


# Read in celltypes TSV
submitter_df <- readr::read_tsv(
    opt$submitter_cell_types_file,
    # read in all columns as character
    col_types = list(.default = readr::col_character()),
    na = character()
  )

# Check columns before proceeding for faster failing:
# TODO: WHAT IF THEY PROVIDE AN ONTOLOGY COLUMN? WE NEED TO PARSE THAT TOO. WILL NEED MORE IF BELOW.
if (!all(c("cell_barcode", "cell_type_assignment") %in% names(submitter_df))) {
  stop("The submitter TSV file must contain columns `cell_barcode` and `cell_type_assignment`.")
}

# Now that we are confident to proceed, read in the sce
sce <- readr::read_rds(opt$sce_file)


# Create submitter_celltype_annotation column
coldata_df <- submitter_df |>
  # filter to relevant library
  dplyr::filter(scpca_library_id == opt$library_id) |>
  # keep columns of interest
  # TODO: WE NEED TO ALLOW FOR AN ONTOLOGY COLUMN
  dplyr::select(
    barcodes = cell_barcode,
    submitter_celltype_annotation = cell_type_assignment
  ) |>
  # join with colData
  dplyr::right_join(
    as.data.frame(colData(sce)),
    by = "barcodes"
  ) |>
  # make any NA values induced by joining into "submitter-excluded"
  dplyr::mutate(
    # use dplyr::if_else, not base, to ensure we end up with character only
    submitter_celltype_annotation = dplyr::if_else(
      is.na(submitter_celltype_annotation),
      "Submitter-excluded",
      submitter_celltype_annotation
    )
  ) |>
  dplyr::distinct()

# Check number of rows before sending back into the SCE object
if (!identical(coldata_df$barcodes, sce$barcodes)) {
  stop("Could not add submitter annotations to SCE object. There should only be one annotation per cell.")
}

# Rejoin with colData, making sure we keep rownames
colData(sce) <- DataFrame(
  coldata_df,
  row.names = colData(sce)$barcodes
)

# Indicate that we have submitter celltypes in metadata,
#  saving in same spot as for actual celltyping
metadata(sce)$celltype_methods <- c(metadata(sce)$celltype_methods, "submitter")

# Write SCE back to file
readr::write_rds(sce, opt$sce_file, compress = "gz")
