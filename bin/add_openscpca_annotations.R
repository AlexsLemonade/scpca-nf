#!/usr/bin/env Rscript

# This script adds openscpca annotations, if provided, to an SCE object's colData.

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
    opt_str = c("--openscpca_cell_types_file"),
    type = "character",
    help = "path to json file containing openscpca cell type annotations"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

# file checks
stopifnot(
  "SingleCellExperiment file name must end in .rds" = stringr::str_ends(opt$sce_file, ".rds"),
  "openscpca cell type annotations file not found" = file.exists(opt$openscpca_cell_types_file),
  "openscpca cell type annotations file must be a json" = stringr::str_ends(opt$openscpca_cell_types_file, ".json")
)


# Read in celltypes json and make sure arrays are not nested lists
json_content <- jsonlite::read_json(
  opt$openscpca_cell_types_file,
  simplifyVector = TRUE
)

celltype_cols <- c("barcodes", "openscpca_celltype_annotation", "openscpca_celltype_ontology")
stopifnot(
  "openscpca cell types annotation file does not contain expected variables" = all(celltype_cols %in% names(json_content))
)

# read in the SCE
sce <- readr::read_rds(opt$sce_file)

# Grab celltypes from json as a data frame to join in with the column data
celltypes_df <- json_content |>
  keep_at(celltype_cols) |>
  as.data.frame()


# join with colData.
# noting by using `left_join()` we preserve the correct order
colData(sce) <- colData(sce) |>
  as.data.frame() |>
  dplyr::left_join(
    celltypes_df,
    by = "barcodes"
  ) |>
  # make any NA values induced by joining with openscpca annotated cells
  dplyr::mutate(
    openscpca_celltype_annotation = dplyr::if_else(
      is.na(openscpca_celltype_annotation),
      "openscpca-excluded",
      openscpca_celltype_annotation
    )
  ) |>
  # Rejoin with colData, making sure we keep rownames
  DataFrame(
    row.names = colnames(sce)
  )

# Indicate that we have openscpca celltypes in metadata,
#  saving in same spot as for actual celltyping
metadata(sce)$celltype_methods <- c(metadata(sce)$celltype_methods, "openscpca")

# add openscpca annotation info as a list to metadata
metadata(sce)$openscpca_annotations_info <- list(
  module_name = json_content$module_name,
  openscpca_nf_version = json_content$openscpca_nf_version,
  release_date = json_content$release_date
)

# Write SCE back to file
readr::write_rds(sce, opt$sce_file, compress = "bz2")
