#!/usr/bin/env Rscript

# This script is used to read in and assign these cell types in the annotated RDS file:
# - CellAssign
# - SingleR
# - consensus cell types
# This script additionally counts the number of normal reference cells in the library and
# adds information to the SCE about reference cells

# import libraries
suppressPackageStartupMessages({
  library(optparse)
  library(SingleCellExperiment)
})

# Functions --------------------------------------------------------------------

#' Get reference information from a reference file
#'
#' @param ref_filename reference file name to parse
#' @param extension file extension to consider during parsing
#'
#' @return list of reference file components to include in SCE metadata
get_ref_info <- function(ref_filename, extension, ref_type) {
  ref_info <- ref_filename |>
    basename() |>
    # select everything before the extension
    stringr::word(1, sep = extension) |>
    # create a vector with name, source, version
    stringr::str_split(pattern = "_") |>
    unlist()

  # account for gene_set_version only being in SingleR refs
  if (ref_type == "SingleR") {
    # if the length is 3, we are using an older version of the
    # SingleR reference file that does not have the gene set version or date
    if (length(ref_info) == 3) {
      # add NA for gene set version and date
      ref_info <- c(ref_info, NA, NA)
    }
    names(ref_info) <- c(
      "ref_name",
      "ref_source",
      "ref_version",
      "gene_set_version",
      "date"
    )
  } else if (ref_type == "CellAssign") {
    names(ref_info) <- c("ref_name", "ref_source", "ref_version")
  }

  return(ref_info)
}

#' Add SingleR results to SCE object
#'
#' @param sce SingleCellExperiment object
#' @param singler_results SingleR results object
#' @param singler_model_file Path to SingleR model file
#'
#' @return Updated SCE object with SingleR annotations
add_singler_results <- function(sce, singler_results, singler_model_file) {
  # get label type from metadata of singler object
  label_type <- metadata(singler_results)$reference_label

  # create a tibble with annotations and barcode
  # later we'll add the annotations into colData by joining on barcodes column
  annotations_df <- tibble::tibble(
    barcodes = rownames(singler_results),
    singler_celltype_annotation = singler_results$pruned.labels,
  )

  # map ontology labels to cell type names, as needed
  # label type will be label.ont if ontology is present
  if (label_type == "label.ont") {
    # end up with columns: barcode, singler_celltype_annotation, singler_celltype_ontology
    annotations_df <- annotations_df |>
      dplyr::left_join(
        # column names: ontology_id, ontology_cell_names
        metadata(singler_results)$cell_ontology_df,
        by = c("singler_celltype_annotation" = "ontology_id")
      ) |>
      # rename columns
      dplyr::rename(
        # ontology should contain the original pruned labels
        singler_celltype_ontology = singler_celltype_annotation,
        # annotation contains the cell names associated with the ontology
        singler_celltype_annotation = ontology_cell_names
      )
  }

  # account for any unclassified cells
  # we need to do this before joining with SingleR results, because otherwise
  # we won't be able to distinguish NA from SingleR vs. NA because it isn't present in SingleR results
  # get list of barcodes in SCE but not in SingleR results
  missing_barcodes <- setdiff(colnames(sce), annotations_df$barcodes)

  # only if there are missing barcodes, append them to the annotations df
  if (length(missing_barcodes) > 0) {
    # create a data frame with unclassified barcodes
    unclassified_df <- data.frame(
      barcodes = missing_barcodes,
      singler_celltype_annotation = "Unclassified cell"
    )

    # combine into one data frame with classified and unclassified cells
    annotations_df <- dplyr::bind_rows(annotations_df, unclassified_df)
  }

  # add annotations to colData
  new_coldata <- colData(sce) |>
    as.data.frame() |>
    dplyr::left_join(annotations_df, by = c("barcodes")) |>
    DataFrame(row.names = colData(sce)$barcodes)
  colData(sce) <- new_coldata

  # get reference name, source and version
  singler_ref_info <- get_ref_info(
    ref_filename = singler_model_file,
    extension = "_model\\.rds",
    ref_type = "SingleR"
  )

  # add singler info to metadata
  metadata(sce)$singler_results <- singler_results
  metadata(sce)$singler_reference <- singler_ref_info[["ref_name"]]
  metadata(sce)$singler_reference_label <- label_type
  metadata(sce)$singler_reference_source <- singler_ref_info[["ref_source"]]
  metadata(sce)$singler_reference_version <- singler_ref_info[["ref_version"]]
  metadata(sce)$singler_gene_set_version <- singler_ref_info[[
    "gene_set_version"
  ]]
  metadata(sce)$singler_date <- singler_ref_info[["date"]]

  return(sce)
}

#' Add CellAssign results to SCE object
#'
#' @param sce SingleCellExperiment object
#' @param predictions CellAssign predictions data frame
#' @param cellassign_ref_file Path to CellAssign reference file
#' @param celltype_ref_metafile Path to cell type reference metadata file
#' @param panglao_ref_file Path to Panglao ontology reference file
#'
#' @return Updated SCE object with CellAssign annotations
add_cellassign_results <- function(
  sce,
  predictions,
  cellassign_ref_file,
  celltype_ref_metafile,
  panglao_ref_file
) {
  # read in panglao ontology reference
  panglao_ref_df <- readr::read_tsv(panglao_ref_file)
  # get cell type with maximum prediction value for each cell
  celltype_assignments <- predictions |>
    # account for the fact that we could have lost some cells we had previously when re-processing through filtering steps
    dplyr::filter(barcode %in% colnames(sce)) |>
    tidyr::pivot_longer(
      !barcode,
      names_to = "celltype",
      values_to = "prediction"
    ) |>
    dplyr::group_by(barcode) |>
    dplyr::slice_max(prediction, n = 1) |>
    dplyr::ungroup() |>
    # add ontology ID column
    dplyr::left_join(panglao_ref_df, by = c("celltype" = "panglao_cell_type"))

  # join by barcode to make sure assignments are in the right order
  celltype_assignments <- data.frame(barcode = sce$barcodes) |>
    dplyr::left_join(celltype_assignments, by = "barcode") |>
    # any cells that are NA were not classified by cellassign
    dplyr::mutate(
      celltype = ifelse(is.na(celltype), "Unclassified cell", celltype)
    )

  # add cell type and prediction to colData
  sce$cellassign_celltype_annotation <- celltype_assignments$celltype
  sce$cellassign_celltype_ontology <- celltype_assignments$ontology
  sce$cellassign_max_prediction <- celltype_assignments$prediction

  # get reference name, source and version
  cellassign_ref_info <- get_ref_info(
    ref_filename = cellassign_ref_file,
    extension = "\\.tsv",
    ref_type = "CellAssign"
  )

  # add entire predictions matrix, ref name, and version to metadata
  metadata(sce)$cellassign_predictions <- predictions
  metadata(sce)$cellassign_reference <- cellassign_ref_info[["ref_name"]]
  metadata(sce)$cellassign_reference_source <- cellassign_ref_info[[
    "ref_source"
  ]]
  metadata(sce)$cellassign_reference_version <- cellassign_ref_info[[
    "ref_version"
  ]]

  # add cellassign reference organs to metadata
  cellassign_organs <- celltype_ref_metafile |>
    readr::read_tsv() |>
    dplyr::filter(celltype_ref_name == cellassign_ref_info[["ref_name"]]) |>
    dplyr::pull(organs)

  if (cellassign_organs == "" | is.na(cellassign_organs)) {
    stop("Failed to obtain CellAssign reference organ list.")
  }
  metadata(sce)$cellassign_reference_organs <- cellassign_organs

  return(sce)
}

#' Add SCimilarity results to SCE object
#'
#' @param sce SingleCellExperiment object
#' @param scimilarity_results_file Path to SCimilarity results file
#' @param scimilarity_model_dir Path to SCimilarity model directory
#'
#' @return Updated SCE object with SCimilarity annotations
add_scimilarity_results <- function(
  sce,
  scimilarity_results_file,
  scimilarity_model_dir
) {
  # read in cell types
  celltype_assignments <- readr::read_tsv(scimilarity_results_file) |>
    # account for the fact that we could have lost some cells we had previously when re-processing through filtering steps
    dplyr::filter(barcode %in% colnames(sce)) |>
    # select the columns to include in the processed object and make sure they are named correctly
    dplyr::select(
      barcodes = barcode,
      scimilarity_celltype_annotation,
      scimilarity_celltype_ontology,
      scimilarity_min_distance = min_dist
    )

  # join by barcode to make sure assignments are in the right order
  celltype_df <- data.frame(barcodes = sce$barcodes) |>
    dplyr::left_join(celltype_assignments, by = "barcodes") |>
    # any cells that are NA were not classified by scimilarity
    dplyr::mutate(
      scimilarity_celltype_annotation = ifelse(
        is.na(scimilarity_celltype_annotation),
        "Unclassified cell",
        scimilarity_celltype_annotation
      )
    )

  # add cell types to colData
  sce$scimilarity_celltype_annotation <- celltype_df$scimilarity_celltype_annotation
  sce$scimilarity_celltype_ontology <- celltype_df$scimilarity_celltype_ontology
  sce$scimilarity_min_distance <- celltype_df$scimilarity_min_distance

  # add model name to metadata
  metadata(sce)$scimilarity_model <- scimilarity_model_dir

  return(sce)
}

#' Assign consensus cell types to SCE object
#'
#' @param sce SingleCellExperiment object
#' @param consensus_celltype_ref Path to consensus cell type reference file
#' @param automated_methods Vector of automated methods used
#'
#' @return Updated SCE object with consensus cell type annotations
assign_consensus_celltypes <- function(
  sce,
  consensus_celltype_ref,
  automated_methods
) {
  # column map indicating which consensus column should be used with which combination of methods
  ref_column_map <- list(
    "consensus" = c("singler", "cellassign", "scimilarity"),
    "cellassign_singler_pair" = c("singler", "cellassign"),
    "singler_scimilarity_pair" = c("singler", "scimilarity"),
    "cellassign_scimilarity_pair" = c("cellassign", "scimilarity")
  )

  # find the appropriate column to use
  ref_column_prefix <- ref_column_map |>
    purrr::keep(\(x) setequal(x, automated_methods)) |>
    names()

  stopifnot(
    "Error getting reference column prefix" = length(ref_column_prefix) == 1
  )

  # define the columns to use in joining based on the existing methods
  join_columns <- glue::glue("{automated_methods}_celltype_ontology")

  consensus_ref_df <- readr::read_tsv(consensus_celltype_ref) |>
    # select columns to use for joining and consensus assignments
    # first make sure the names match what we expect
    dplyr::rename(
      cellassign_celltype_ontology = panglao_ontology,
      singler_celltype_ontology = blueprint_ontology,
      scimilarity_celltype_ontology = scimilarity_ontology
    ) |>
    # now just filter to join columns and ref column with consensus cell types
    dplyr::select(dplyr::all_of(join_columns), dplyr::starts_with(ref_column_prefix)) |>
    # only keep unique combos
    dplyr::distinct() |>
    # make sure the columns used to get the consensus cell type actually have the consensus_ prefix rather than singler_cellassign_pair_, etc.
    dplyr::rename_with(
      \(x) stringr::str_replace(x, ref_column_prefix, "consensus"),
      .cols = dplyr::starts_with(ref_column_prefix)
    )

  # create df with consensus assignments
  celltype_df <- colData(sce) |>
    as.data.frame() |>
    dplyr::select(
      barcodes,
      dplyr::contains("celltype") # get any available cell type columns with ontology
    ) |>
    # then add consensus labels
    dplyr::left_join(
      consensus_ref_df,
      by = join_columns,
      relationship = "many-to-one" # account for multiple of the same cell type
    ) |>
    # use unknown for NA annotation but keep ontology ID as NA
    tidyr::replace_na(list(consensus_annotation = "Unknown"))

  # add consensus cell type and ontology to sce
  sce$consensus_celltype_annotation <- celltype_df$consensus_annotation
  sce$consensus_celltype_ontology <- celltype_df$consensus_ontology

  # add a note to metadata indicating which methods were used to assign consensus
  metadata(sce)$consensus_celltype_methods <- automated_methods

  return(sce)
}


#' Assign infercnv_status to SCE metadata
#'
#' If the final infercnv_status is "", no edge cases were encountered
#' If it can be determined, this function will also add infercnv_diagnosis_groups to the SCE metadata
#'
#' @param sce SingleCellExperiment object
#' @param diagnosis_celltype_ref Path to tsv mapping broad diagnoses to consensus validation groups
#' @param diagnosis_groups_ref Path to tsv mapping broad diagnoses to individual diagnoses
#'
#' @returns Updated SCE object with infercnv_status and infercnv_diagnosis_groups added to SCE metadata
assign_infercnv_status <- function(
  sce,
  diagnosis_groups_ref,
  diagnosis_celltype_ref
) {
  # unique in case of multiplexed
  diagnosis <- unique(metadata(sce)$sample_metadata$diagnosis)

  # Assign status and return early if there is no usable diagnosis information
  # these edge cases will not have a broad_diagnosis saved in the SCE metadata
  if (is.null(diagnosis)) {
    metadata(sce)$infercnv_status <- "unknown_diagnosis"
    return(sce)
  }
  if (length(diagnosis) > 1 && !file.exists(diagnosis_groups_ref)) {
    metadata(sce)$infercnv_status <- "multiple_diagnoses_multiplexed"
    return(sce)
  }

  # Find the broad_diagnosis
  if (!file.exists(diagnosis_groups_ref)) {
    broad_diagnosis <- diagnosis # if no map file, use the given diagnosis
  } else {
    diagnosis_groups <- readr::read_tsv(diagnosis_groups_ref)
    broad_diagnosis <- data.frame(sample_diagnosis = diagnosis) |>
      dplyr::left_join(diagnosis_groups, by = "sample_diagnosis") |>
      # replace NA diagnosis_group with sample diagnosis
      # this means diagnosis_group will never be length 0
      dplyr::mutate(
        diagnosis_group = dplyr::coalesce(diagnosis_group, sample_diagnosis)
      ) |>
      dplyr::pull(diagnosis_group) |>
      # in case of multiplexed
      unique()
  }
  # update broad_diagnosis - remove non-cancerous for multiplexed edge case
  if (length(broad_diagnosis) > 1 & "Non-cancerous" %in% broad_diagnosis) {
    broad_diagnosis <- broad_diagnosis[!broad_diagnosis == "Non-cancerous"]
  }

  # save broad diagnosis to SCE before remaining checks
  metadata(sce)$infercnv_diagnosis_groups <- broad_diagnosis


  # Assign remaining statuses for edge cases, returning early to reduce nesting
  if (length(broad_diagnosis) > 1) {
    metadata(sce)$infercnv_status <- "multiple_diagnosis_groups_multiplexed"
    return(sce)
  }
  if (broad_diagnosis == "Non-cancerous") { # at this point we know it's length 1
    metadata(sce)$infercnv_status <- "skipped_non_cancerous"
    return(sce)
  }
  if (!file.exists(diagnosis_celltype_ref)) {
    metadata(sce)$infercnv_status <- "no_diagnosis_celltype_reference"
    return(sce)
  } else {
    diagnosis_celltype_df <- readr::read_tsv(diagnosis_celltype_ref)
    if (!(broad_diagnosis %in% diagnosis_celltype_df$diagnosis_group)) {
      metadata(sce)$infercnv_status <- "unknown_reference_celltypes"
      return(sce)
    }
  }

  # If we made it here, no edge case was found
  metadata(sce)$infercnv_status <- ""
  return(sce)
}


#' Add reference cell information for inferCNV to an SCE
#'
#'
#' @param sce SingleCellExperiment object
#' @param consensus_validation_ref Path to tsv mapping consensus validation groups to consensus labels
#' @param diagnosis_celltype_ref Path to tsv mapping broad diagnoses to consensus validation groups
#'
#' @returns Updated SCE object with reference cell count information
add_infercnv_reference_cells <- function(
  sce,
  consensus_validation_ref,
  diagnosis_celltype_ref
) {
  consensus_validation_df <- readr::read_tsv(opt$consensus_validation_ref)
  diagnosis_celltype_df <- readr::read_tsv(opt$diagnosis_celltype_ref)

  # get the cell type groups to consider for this diagnosis
  reference_validation_groups <- diagnosis_celltype_df |>
    dplyr::filter(diagnosis_group == metadata(sce)$infercnv_diagnosis_groups) |>
    tidyr::separate_longer_delim(celltype_groups, delim = ",") |>
    dplyr::pull(celltype_groups) |>
    # remove any leading or trailing spaces
    stringr::str_trim()

  # get the consensus cell types
  ref_df <- consensus_validation_df |>
    dplyr::filter(validation_group_annotation %in% reference_validation_groups) |>
    dplyr::select(consensus_ontology, consensus_annotation) |>
    dplyr::distinct()

  # Add reference cell information to SCE
  metadata(sce)$infercnv_reference_celltypes <- ref_df$consensus_annotation # vector of reference cell types
  sce$is_infercnv_reference <- sce$consensus_celltype_ontology %in% ref_df$consensus_ontology

  # include the total number of reference cells in the metadata
  metadata(sce)$infercnv_num_reference_cells <- sum(sce$is_infercnv_reference)

  return(sce)
}


# Main script ----------------------------------------------------------------
option_list <- list(
  make_option(
    opt_str = c("-i", "--input_sce_file"),
    type = "character",
    help = "path to rds file with input sce object",
    default = ""
  ),
  make_option(
    opt_str = c("-o", "--output_sce_file"),
    type = "character",
    help = "path to output rds file to store annotated sce object. Must end in .rds",
    default = ""
  ),
  make_option(
    opt_str = c("--singler_results"),
    type = "character",
    help = "path to rds file containing SingleR results object",
    default = ""
  ),
  make_option(
    opt_str = c("--singler_model_file"),
    type = "character",
    help = "Name of file containing a single model generated for SingleR annotation.
            File name is expected to be in form: `<ref_name>_<source>_<version>_<gene_set_version>_<date>_model.rds`.",
    default = ""
  ),
  make_option(
    opt_str = c("--cellassign_predictions"),
    type = "character",
    help = "path to tsv file containing the prediction matrix returned by running CellAssign",
    default = ""
  ),
  make_option(
    opt_str = c("--cellassign_ref_file"),
    type = "character",
    help = "Name of marker by cell type reference file used with CellAssign.
            File name is expected to be in form: `<ref_name>_<source>_<version>.tsv`",
    default = ""
  ),
  make_option(
    opt_str = c("--celltype_ref_metafile"),
    type = "character",
    help = "Metadata TSV file containing cell type reference metadata.
      This file is used to obtain a list of organs used for CellAssign annotation and is only required if CellAssign results provided as input.",
    default = ""
  ),
  make_option(
    opt_str = c("--panglao_ontology_ref"),
    type = "character",
    help = "Path to TSV file with panglao assignments and associated cell ontology ids.
      This file is used to assign Cell Ontology identifier to the CellAssign annotations",
    default = ""
  ),
  make_option(
    opt_str = c("--scimilarity_results"),
    type = "character",
    help = "path to tsv file containing SCimilarity results",
    default = ""
  ),
  make_option(
    opt_str = c("--scimilarity_model_dir"),
    type = "character",
    help = "Name of directory with the model used for SCimilarity",
    default = ""
  ),
  make_option(
    opt_str = c("--consensus_celltype_ref"),
    type = "character",
    help = "Path to file containing the reference for assigning consensus cell type labels",
    default = ""
  ),
  make_option(
    opt_str = c("--consensus_validation_ref"),
    type = "character",
    help = "Path to TSV file mapping consensus validation groups to consensus labels for counting normal reference cells intended as input to inferCNV",
    default = ""
  ),
  make_option(
    opt_str = c("--diagnosis_celltype_ref"),
    type = "character",
    help = "Path to TSV file mapping broad diagnoses to consensus validation groups for counting normal reference cells intended as input to inferCNV.
      This file should have columns `diagnosis_group` and `celltype_groups` where the latter is a comma-separated list of consensus validation groups",
    default = ""
  ),
  make_option(
    opt_str = c("--diagnosis_groups_ref"),
    type = "character",
    help = "Path to TSV file mapping broad diagnoses to individual diagnoses for counting normal reference cells intended as input to inferCNV.
      This file should have columns `diagnosis_group` for the broad diagnosis and `sample_diagnosis` with individual diagnoses in ScPCA",
    default = ""
  ),
  make_option(
    opt_str = c("--reference_cell_count_file"),
    type = "character",
    help = "Path to write number of calculated inferCNV reference cells to.
      This calculation is only performed if `diagnosis_celltype_ref`, `diagnosis_groups_ref`, and `consensus_validation_ref` are provided.
      If not calculated, this file will be empty",
    default = "reference_cell_count.txt"
  ),
  make_option(
    opt_str = c("--reference_cell_hash_file"),
    type = "character",
    help = "Path to write out unique hash for all concatenated barcodes in reference cell set; used for checkpointing.
      If not calculated, this file will be empty",
    default = "reference_cell_hash.txt"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

# check that input file exists, output file ends in rds, and count files are provided
stopifnot(
  "Missing input SCE file" = file.exists(opt$input_sce_file),
  "output sce file name must end in .rds" = stringr::str_ends(opt$output_sce_file, ".rds"),
  "output file to store counted normal cells was not provided" =
    !is.null(opt$reference_cell_count_file) && opt$reference_cell_count_file != "",
  "output file to store reference cell hash was not provided" =
    !is.null(opt$reference_cell_hash_file) && opt$reference_cell_hash_file != ""
)

# read in input files
sce <- readr::read_rds(opt$input_sce_file)

# We'll count inferCNV reference cells when assigning consensus cell types
reference_cell_count <- ""
reference_cell_hash <- ""

# set automated methods list
# we'll use this to assign consensus and add cell type methods to the metadata
automated_methods <- c()

# SingleR results --------------------------------------------------------------

has_singler <- file.exists(opt$singler_results)
if (has_singler) {
  # check singler model has been provided
  stopifnot(
    "Singler model filename must be provided" = opt$singler_model_file != ""
  )

  singler_results <- readr::read_rds(opt$singler_results)

  # add singler results to sce
  sce <- add_singler_results(sce, singler_results, opt$singler_model_file)

  # note that singler is included in methods
  automated_methods <- c(automated_methods, "singler")
}

# CellAssign results -----------------------------------------------------------

has_cellassign <- file.exists(opt$cellassign_predictions)
if (has_cellassign) {
  # check that cellassign reference info is provided
  stopifnot(
    "CellAssign reference filename must be provided" = opt$cellassign_ref_file != "",
    "Cell type reference metadata file does not exist" = file.exists(opt$celltype_ref_metafile),
    "Panglo ontology reference file does not exist" = file.exists(opt$panglao_ontology_ref)
  )

  # if cell assign predictions file exists but is empty then cell assign failed and we want to account for that with Not Run
  if (file.size(opt$cellassign_predictions) > 0) {
    # read in predictions file
    cellassign_df <- readr::read_tsv(opt$cellassign_predictions)
  } else {
    cellassign_df <- NULL
    has_cellassign <- FALSE # reset to false so that we don't add in consensus cell types
  }

  # if the only column is the barcode column or if the predictions file was empty
  # then CellAssign didn't complete successfully
  # otherwise add in cell type annotations and metadata to SCE
  if (is.null(cellassign_df) || all(colnames(cellassign_df) == "barcode")) {
    # if failed then note that in the cell type column
    sce$cellassign_celltype_annotation <- "Not run"
    has_cellassign <- FALSE # reset to false if cellassign didn't complete so we don't add consensus
  } else {
    # add cellassign results to sce
    sce <- add_cellassign_results(
      sce,
      cellassign_df,
      opt$cellassign_ref_file,
      opt$celltype_ref_metafile,
      opt$panglao_ontology_ref
    )

    # add cellassign as celltype method
    automated_methods <- c(automated_methods, "cellassign")
  }
}

# SCimilarity ------------------------------------------------------------------

has_scimilarity <- file.exists(opt$scimilarity_results)
if (has_scimilarity) {
  # check that scimilarity model info is provided
  stopifnot(
    "SCimilarity model directory name must be provided" = opt$scimilarity_model_dir != ""
  )

  # add scimilarity results to sce
  sce <- add_scimilarity_results(
    sce,
    opt$scimilarity_results,
    opt$scimilarity_model_dir
  )

  # add scimilarity as celltype method
  automated_methods <- c(automated_methods, "scimilarity")
}

# add automated methods to the metadata of the object
# note that if `metadata(sce)$celltype_methods` doesn't exist yet, this will
#  come out to just the automated methods
metadata(sce)$celltype_methods <- c(
  metadata(sce)$celltype_methods,
  automated_methods
)

# Consensus assignment ---------------------------------------------------------

# set columns to use for joining based on methods that are present
# define the prefix of the column in the consensus reference that contains the appropriate consensus term given the provided methods
# e.g., all three methods use the main consensus_annotation column
# if the library only has scimilarity and singler, use singler_scimilarity_pair_annotation column, etc.


# if there's at least two methods then assign consensus using the available methods
if (length(automated_methods) > 1) {
  # now make sure that reference file exists
  stopifnot(
    "Consensus cell type reference file does not exist" = file.exists(opt$consensus_celltype_ref)
  )

  # assign consensus cell types
  sce <- assign_consensus_celltypes(
    sce,
    opt$consensus_celltype_ref,
    automated_methods
  )

  # Count the number of reference cells -------------------------
  stopifnot(
    "Consensus cell type validation file does not exist" = file.exists(opt$consensus_validation_ref)
  )

  # assign status, except for case where there are no consensus annotations
  sce <- assign_infercnv_status(
    sce,
    opt$diagnosis_groups_ref,
    opt$diagnosis_celltype_ref
  )
  # calculate if passing status
  if (metadata(sce)$infercnv_status == "") {
    sce <- add_infercnv_reference_cells(
      sce,
      opt$consensus_validation_ref,
      opt$diagnosis_celltype_ref
    )
    reference_cell_count <- metadata(sce)$infercnv_num_reference_cells
    reference_cell_hash <- sort(sce$barcodes[sce$is_infercnv_reference]) |>
      paste(collapse = "") |>
      digest::digest()
  }
} else {
  metadata(sce)$infercnv_status <- "no_consensus"
}

# export annotated object with cell type assignments
readr::write_rds(sce, opt$output_sce_file, compress = "bz2")

# export normal cell count and hash
readr::write_file(reference_cell_count, opt$reference_cell_count_file)
readr::write_file(reference_cell_hash, opt$reference_cell_hash_file)
