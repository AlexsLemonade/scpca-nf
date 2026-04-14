#!/usr/bin/env Rscript

# check_sce_format.R
# Usage: Rscript check_sce_format.R <sce_file> <reference_json> <object_type>
# object_type: one of "unfiltered", "filtered", "processed"

# output is either an empty text file (no errors) or text file with list of errors

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(jsonlite)
  library(optparse)
})

option_list <- list(
  make_option(
    opt_str = c("--sce_file"),
    type = "character",
    default = "",
    help = "Path to RDS file containing a SingleCellExperiment object from scpca-nf"
  ),
  make_option(
    opt_str = c("--object_type"),
    type = "character",
    help = "Type of object from scpca-nf, either unfiltered, filtered, or processed"
  ),
  make_option(
    opt_str = c("--reference_file"),
    type = "character",
    default = "",
    help = "Path to JSON file with formatting reference"
  ),
  make_option(
    opt_str = c("--output_file"),
    type = "character",
    default = "formatting_errors.txt",
    help = ".txt file to save any identified errors"
  )
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

# TODO: Remove testing defaults before merging
opt$sce_file <- "~/Downloads/SCPCL000001_processed.rds"
opt$object_type <- "processed"
opt$reference_file <- "../scpca-nf/references/sce-formatting-reference.json"

# Set up -----------------------------------------------------------------------

# make sure inputs are correct/exist
stopifnot(
  "sce file does not exist" = file.exists(opt$sce_file),
  "Object type must be one of unfiltered, filtered, or processed" = opt$object_type %in% c("unfiltered", "filtered", "processed"),
  "Reference file does not exist" = file.exists(opt$reference_file),
  "output file must end in `.txt`" = stringr::str_ends(opt$gene_exp_output_file, "\\.txt")
)

# read in sce
sce <- readr::read_rds(opt$sce_file)

# read in formatting file and subset to object type
ref_list <- jsonlite::read_json(opt$reference_file)[[opt$object_type]]

# set empty error messages
errors <- ""

# Check assays and ensure counts are rounded -----------------------------------
missing_assays <- assayNames(sce)[!ref_list$assayNames %in% assayNames(sce)]
if (length(missing_assays) > 0) {
  errors <- glue::glue("{errors}
                       Missing assay: {missing_assays}
                       ")
}

all_values <- as.vector(counts(sce))
if (any(all_values != floor(all_values))) {
  errors <- glue::glue("
                       {errors}
                       'counts' assay does not contain rounded integers
                       ")
}

# check colData names and type -------------------------------------------------

check_names_and_types <- function(
  data,
  ref
) {
  match_df <- ref |>
    purrr::imap(\(type, name){
      is_present <- name %in% names(data)
      if (is_present) {
        obs_type <- class(data[[name]])
        is_type_match <- obs_type == type
      } else {
        obs_type <- NULL
      }

      return(
        list(
          "is_present" = is_present,
          "expected_type" = type,
          "observed_type" = obs_type,
          "is_type_match" = is_type_match
        )
      )
    }) |>
    dplyr::bind_rows(.id = "name")

  return(match_df)
}

coldata_match_df <- check_names_and_types(colData(sce), ref_list$colData)
rowdata_match_df <- check_names_and_types(rowData(sce), ref_list$rowData)
metadata_match_df <- check_names_and_types(metadata(sce), ref_list$metadata)

# conditional
conditional_coldata_ref <- ref_list$colData_conditional

# check for existence of specific results
conditionals_vec <- c(
  has_normalization,
  has_infercnv,
  has_clusters = "cluster" %in% names(colData(sce)),
  has_consensus = "consensus_celltype_annotation" %in% names(colData(sce)),
  has_singler = "singler" %in% metadata(sce)$celltype_methods,
  has_cellassign = "cellassign" %in% metadata(sce)$celltype_methods,
  has_scimilarity = "scimilarity" %in% metadata(sce)$celltype_methods,
  has_openscpca = "openscpca" %in% metadata(sce)$celltype_methods,
  has_submitter = "submitter" %in% metadata(sce)$celltype_methods &&
    !all(is.na(sce$submitter_celltype_annotation)) # make sure they aren't all NA
)

conditional_coldata_match_df <- names(conditionals_vec[which(conditionals_vec)]) |>
  purrr::map(\(condition){
    df <- check_names_and_types(colData(sce), conditional_coldata_ref[[condition]])
    df$conditional <- condition
    return(df)
  }) |>
  dplyr::bind_rows()

# check reduced Dim if processed

# TODO: handle altExps

# cat error messages and export file
