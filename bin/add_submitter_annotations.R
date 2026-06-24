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

# Check required columns before proceeding for faster failing:
if (!all(c("cell_barcode", "cell_type_assignment") %in% names(submitter_df))) {
  stop("The submitter TSV file must contain columns `cell_barcode` and `cell_type_assignment`.")
}

# Now that we are confident to proceed, read in the sce
sce <- readr::read_rds(opt$sce_file)

# Identify any extra columns beyond the required ones in the submitters file
# and CL ontology Id which will get its own name if present
# specify any columns that might be present that we don't need to include in the output object (e.g., scpca_sample_id or submitter_id)
extra_cols <- setdiff(
  names(submitter_df),
  c(
    "scpca_library_id",
    "submitter_id",
    "scpca_sample_id",
    "cell_barcode",
    "cell_type_assignment",
    "CL_ontology_id"
  )
)

# Rename extra columns with `submitter_` prefix, replacing dashes/periods/spaces
# with underscores and converting to lower case
renamed_extra_cols <- extra_cols |>
  tolower() |>
  stringr::str_replace_all("[\\-\\.\\s]+", "_") |>
  (\(x) paste0("submitter_", x))()

# Build a named vector for renaming: old name -> new name
extra_col_rename <- setNames(extra_cols, renamed_extra_cols)

# Create submitter_celltype_annotation column
submitter_df <- submitter_df |>
  # filter to relevant library
  dplyr::filter(scpca_library_id == opt$library_id) |>
  # rename specific columns for barcode and cell type assignment
  dplyr::rename(
    barcodes = cell_barcode,
    submitter_celltype_annotation = cell_type_assignment
  ) |>
  # rename extra columns with submitter_ prefix
  dplyr::rename(!!!extra_col_rename) |>
  dplyr::distinct() 

# specifically rename CL ontology id if present
if ("CL_ontology_id" %in% colnames(submitter_df)) {
  submitter_df <- submitter_df |>
    # select only the columns we want to include and rename ontology ID
    dplyr::select(
      barcodes,
      submitter_celltype_annotation,
      submitter_celltype_ontology = CL_ontology_id,
      renamed_extra_cols
    )

  # assign the new column to a variable to input to the list of columns
  # to replace NA with submitter-excluded
  ontology_column <- "submitter_celltype_ontology"
} else {
  submitter_df <- submitter_df |>
    dplyr::select(
      barcodes,
      submitter_celltype_annotation,
      renamed_extra_cols
    )

  ontology_column <- NULL
}

# All submitter columns
submitter_annotation_cols <- c("submitter_celltype_annotation", ontology_column, renamed_extra_cols)

# join with colData.
# noting by using `left_join()` we preserve the correct order
coldata_df <- colData(sce) |>
  as.data.frame() |>
  dplyr::left_join(
    submitter_df,
    by = "barcodes"
  ) |>
  # for cells that are not present in the submitter file
  # fill in values in submitter columns with "Submitter-excluded"
  dplyr::mutate(
    dplyr::across(
      dplyr::all_of(submitter_annotation_cols),
      # use dplyr::if_else, not base, to ensure we end up with character only
      \(x) dplyr::if_else(!barcodes %in% submitter_df$barcodes, "Submitter-excluded", x)
    )
  )


# Check that barcodes are correct before sending back into the SCE object
if (!identical(coldata_df$barcodes, sce$barcodes)) {
  stop("Failed to add submitter annotations to SCE object.")
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
readr::write_rds(sce, opt$sce_file, compress = "bz2")
