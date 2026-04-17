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

# Functions --------------------------------------------------------------------

#' Check that expected names and types are present in a data slot
#'
#' @param data A data slot from an SCE object (e.g., colData, rowData, metadata)
#' @param ref A named list mapping expected slots to their expected contents and types
#'
#' @return A data frame with columns: name, is_present, expected_type,
#'   observed_type, and is_type_match
check_names_and_types <- function(data, ref) {
  ref |>
    purrr::imap(\(type, name) {
      is_present <- name %in% names(data)

      if (!is_present) {
        obs_type <- NA
        is_type_match <- NA
      } else {
        # first get the observed type
        obs_type <- class(data[[name]])

        # now see if it matches the expected type
        if (name == "sample_type" && is(data[[name]], "list")) {
          # sample type may be a list for multiplexed
          is_type_match <- is(data[[name]][[0]], type)
        } else {
          is_type_match <- is(data[[name]], type)
        }
      }

      list(
        "is_present"     = is_present,
        "expected_type"  = type,
        "observed_type"  = obs_type,
        "is_type_match"  = is_type_match
      )
    }) |>
    dplyr::bind_rows(.id = "name")
}


#' Check names and types for conditional data slot entries
#'
#' For each condition that is TRUE for the given object, checks that the
#' corresponding expected names and types are present in the data slot.
#'
#' @param conditions_present A character vector of condition names that are
#'   TRUE for the given object (i.e., the subset of names from conditionals_vec
#'   where the value is TRUE)
#' @param data A data slot from an SCE object (e.g., colData, metadata)
#' @param ref A named list where each element corresponds to a condition name
#'   and contains a named list of expected column names and types
#'
#' @return A data frame with the same columns as check_names_and_types, plus
#'   a `conditional` column indicating which condition triggered the check
check_conditional_names_and_type <- function(conditions_present, data, ref) {
  conditions_present |>
    purrr::map(\(condition) {
      df <- check_names_and_types(data, ref[[condition]])
      df$conditional <- condition
      return(df)
    }) |>
    dplyr::bind_rows()
}


#' Check required and conditional col/row/metadata slots of an SCE object
#'
#' Runs both required and conditional checks against a reference list,
#' then combines the results into a single data frame.
#'
#' @param sce_exp A SingleCellExperiment exp (can be a main or alt experiment)
#' @param ref_list A named list containing the reference for the object type,
#'   with elements: colData, rowData, colData_conditional, metadata_conditional
#' @param conditionals_vec A named logical vector indicating which conditions
#'   are TRUE for the object being checked. Should be defined using the
#'   top-level SCE (not an altExp) since many conditions depend on top-level
#'   metadata (e.g., altExpNames, celltype_methods)
#'
#' @return A data frame combining required and conditional check results,
#'   with a `slot` column indicating the source slot and a `conditional` column
#'   (NA for required checks)
check_sce_object <- function(sce_exp, ref_list, conditionals_vec) {
  ############ required colData and rowData checks ############
  match_df <- list(
    # use square brackets to avoid json doing any partial matching
    "colData" = list(data = colData(sce_exp), ref = ref_list[["colData"]]),
    "rowData" = list(data = rowData(sce_exp), ref = ref_list[["rowData"]]),
    "metadata" = list(data = metadata(sce_exp), ref = ref_list[["metadata"]])
  ) |>
    # skip slots that have no reference defined (e.g. row/colData absent for altExps)
    purrr::keep(\(x) length(x$ref) > 0) |>
    purrr::map(\(x) check_names_and_types(x$data, x$ref)) |>
    dplyr::bind_rows(.id = "slot")

  ####### conditional colData and metadata checks #######
  # only check conditions that are true for this object
  true_conditions <- names(conditionals_vec[which(conditionals_vec)])

  if (length(true_conditions) > 1) {
    conditional_match_df <- list(
      colData  = list(data = colData(sce_exp), ref = ref_list[["colData_conditional"]]),
      metadata = list(data = metadata(sce_exp), ref = ref_list[["metadata_conditional"]])
    ) |>
      # skip slots that have no reference defined
      purrr::keep(\(x) length(x$ref) > 0) |>
      purrr::map(\(x) check_conditional_names_and_type(
        intersect(true_conditions, names(x$ref)),
        x$data,
        x$ref
      )) |>
      dplyr::bind_rows(.id = "slot")

    # combine required and conditional checks
    all_match_df <- match_df |>
      dplyr::mutate(conditional = NA) |>
      dplyr::bind_rows(conditional_match_df)
  } else {
    all_match_df <- match_df
  }

  return(all_match_df)
}


#' Check that all expected assays are present in an SCE object
#'
#' @param sce_exp A SingleCellExperiment exp (can be a main or alt experiment)
#' @param ref_assay_names A character vector of expected assay names
#' @param label A string used in error messages to identify the object
#'   (e.g., "main SCE", "adt altExp"). Default is "SCE".
#'
#' @return A character vector of error messages for any missing assays,
#'   or an empty character vector if all assays are present
check_assays <- function(sce_exp, ref_assay_names, label = "SCE") {
  missing_assays <- ref_assay_names[!ref_assay_names %in% assayNames(sce_exp)]

  if (length(missing_assays) > 0) {
    glue::glue("Missing assay '{missing_assays}' from {label}.")
  } else {
    character(0)
  }
}


#' Check col/row/metadata match results and append errors to the error string
#'
#' Takes the output of check_sce_object and adds formatted error messages
#' for any missing fields or type mismatches.
#'
#' @param match_df A data frame returned by check_sce_object
#' @param errors A string of accumulated error messages to append to
#' @param label A string used in error messages to identify the object
#'   (e.g., "main SCE", "adt altExp"). Default is "SCE".
#'
#' @return An updated error string with any new errors appended
collect_match_errors <- function(match_df, errors, label = "SCE") {
  # TODO: Remove this if we accept using NA_real_
  # remove double prob_compromised_cutoff if its present for conditional and regular metadata
  # this value is present with or without miQC with a different type so we need to keep both checks in the reference
  num_prob_compromised_checks <- match_df$name[stringr::str_detect(match_df$name, "prob_compromised_cutoff")]
  if (length(num_prob_compromised_checks) == 2) {
    match_df <- match_df |>
      dplyr::filter(!(name == "prob_compromised_cutoff" & is.na(conditional)))
  }

  # missing field errors
  missing_df <- dplyr::filter(match_df, !is_present)
  if (nrow(missing_df) > 0) {
    errors <- c(
      errors,
      glue::glue("
      Missing '{missing_df$name}' from {label} {missing_df$slot}.
        ")
    )
  }

  # type mismatch errors
  mismatch_df <- dplyr::filter(match_df, is_present, !is_type_match)
  if (nrow(mismatch_df) > 0) {
    errors <- c(
      errors,
      glue::glue("
      Type mismatch in '{mismatch_df$name}' from {label} {mismatch_df$slot}.
      Expected {mismatch_df$expected_type}, but found {mismatch_df$observed_type}.
      ")
    )
  }

  return(errors)
}


# Arguments --------------------------------------------------------------------

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

opt <- parse_args(OptionParser(option_list = option_list))

# Set up -----------------------------------------------------------------------

stopifnot(
  "sce file does not exist" = file.exists(opt$sce_file),
  "Object type must be one of unfiltered, filtered, or processed" = opt$object_type %in% c("unfiltered", "filtered", "processed"),
  "Reference file does not exist" = file.exists(opt$reference_file),
  "output file must end in `.txt`" = stringr::str_ends(opt$output_file, "\\.txt")
)

sce <- readr::read_rds(opt$sce_file)
all_ref_list <- jsonlite::read_json(opt$reference_file)[[opt$object_type]]

# initialize empty error string
errors <- ""

# Check main SCE assays --------------------------------------------------------

errors <- c(errors, check_assays(sce, all_ref_list$assayNames, label = "main SCE"))

# check that counts assay contains rounded integers
if (all(counts(sce)@x == floor(counts(sce)@x))) {
  errors <- glue::glue("{errors}
                       'counts' assay does not contain rounded integers.
                       ")
}

# Define conditionals for the main SCE object ----------------------------------
# these are evaluated once on the top-level sce and reused for altExp checks

has_adt <- "adt" %in% altExpNames(sce)
has_cellhash <- "cellhash" %in% altExpNames(sce)

conditionals_vec <- c(
  # mapping tools
  "alevin-fry" = metadata(sce)$mapping_tool == "alevin-fry",
  "cellranger-multi" = metadata(sce)$mapping_tool == "cellranger multi",

  # preprocessing
  umi_filtering = metadata(sce)$filtering_method == "UMI cutoff",
  # the only way to confirm miQC was run successfully is if this cutoff is present and numeric
  # TODO: use the has_miQC entry instead
  has_miQC = class(metadata(sce)$prob_compromised_cutoff) == "numeric",
  has_normalization = metadata(sce)$normalization == "normalization",

  # clustering and cell typing
  has_clusters = "cluster" %in% names(colData(sce)),
  has_consensus = "consensus_celltype_annotation" %in% names(colData(sce)),
  has_singler = "singler" %in% metadata(sce)$celltype_methods,
  has_cellassign = "cellassign" %in% metadata(sce)$celltype_methods,
  has_scimilarity = "scimilarity" %in% metadata(sce)$celltype_methods,
  has_openscpca = "openscpca" %in% metadata(sce)$celltype_methods,
  # submitter annotations are only valid if not all NA
  has_submitter = "submitter" %in% metadata(sce)$celltype_methods &&
    !all(is.na(sce$submitter_celltype_annotation)),
  has_infercnv = !is.null(metadata(sce)$infercnv_status),

  # additional modalities
  has_adt = has_adt,
  has_cellhash = has_cellhash,
  has_hashedDrops =
    any(stringr::str_detect(colnames(colData(sce)), "hashedDrops")),
  has_HTODemux =
    any(stringr::str_detect(colnames(colData(sce)), "HTODemux")),
  has_vireo =
    any(stringr::str_detect(colnames(colData(sce)), "vireo"))
)

# Check main SCE col/row/metadata ----------------------------------------------

all_match_df <- check_sce_object(sce, all_ref_list, conditionals_vec)
errors <- collect_match_errors(all_match_df, errors, label = "main SCE")

# Check ADT altExp -------------------------------------------------------------

if (has_adt) {
  adt_sce <- altExp(sce, "adt")
  adt_ref <- all_ref_list$altExp$adt

  adt_conditionals_vec <- c(
    "has_negative_control" = any(rowData(adt_sce)$target_type %in% c("neg_control")),
    "no_negative_control" = all(rowData(adt_sce)$target_type %in% c("pos_control", "target")),
    "alevin-fry" = metadata(sce)$mapping_tool == "alevin-fry",
    "cellranger-multi" = metadata(sce)$mapping_tool == "cellranger multi"
  )

  adt_match_df <- check_sce_object(adt_sce, adt_ref, adt_conditionals_vec)

  # check assay errors first
  errors <- c(errors, check_assays(adt_sce, adt_ref$assayNames, label = "adt altExp"))
  # now add in col/row/metadata
  errors <- collect_match_errors(adt_match_df, errors, label = "adt altExp")
}

# Check cellhash altExp --------------------------------------------------------

if (has_cellhash) {
  cellhash_sce <- altExp(sce, "cellhash")
  cellhash_ref <- all_ref_list$altExp$cellhash

  cellhash_conditionals_vec <- c(
    has_hashedDrops =
      any(stringr::str_detect(colnames(colData(cellhash_sce)), "hashedDrops")),
    has_HTODemux =
      any(stringr::str_detect(colnames(colData(cellhash_sce)), "HTODemux"))
  )

  cellhash_match_df <- check_sce_object(cellhash_sce, cellhash_ref, cellhash_conditionals_vec)

  # check assay errors first
  errors <- c(errors, check_assays(cellhash_sce, cellhash_ref$assayNames, label = "cellhash altExp"))
  # now check col/row/metadata
  errors <- collect_match_errors(cellhash_match_df, errors, label = "cellhash altExp")
}

# Check reducedDims (processed objects only) -----------------------------------

if (opt$object_type == "processed") {
  if (all(reducedDimNames(sce) != all_ref_list$reducedDimNames)) {
    errors <- glue::glue(
      "{errors}
       reducedDimNames do not match expected: {all_ref_list$reducedDimNames}.
      "
    )
  }
}

# Write output -----------------------------------------------------------------

writeLines(errors, opt$output_file)
