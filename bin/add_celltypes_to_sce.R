#!/usr/bin/env Rscript

# This script is used to read in the predictions from CellAssign
# and assign cell types in the annotated RDS file

# import libraries
suppressPackageStartupMessages({
  library(optparse)
  library(SingleCellExperiment)
})

# set up arguments
option_list <- list(
  make_option(
    opt_str = c("-i", "--input_sce_file"),
    type = "character",
    help = "path to rds file with input sce object"
  ),
  make_option(
    opt_str = c("-o", "--output_sce_file"),
    type = "character",
    help = "path to output rds file to store annotated sce object. Must end in .rds"
  ),
  make_option(
    opt_str = c("--singler_results"),
    type = "character",
    help = "path to rds file containing SingleR results object"
  ),
  make_option(
    opt_str = c("--singler_model_file"),
    type = "character",
    help = "Name of file containing a single model generated for SingleR annotation.
            File name is expected to be in form: `<ref_name>_<source>_<version>_model.rds`."
  ),
  make_option(
    opt_str = c("--cellassign_predictions"),
    type = "character",
    help = "path to tsv file containing the prediction matrix returned by running CellAssign"
  ),
  make_option(
    opt_str = c("--cellassign_ref_file"),
    type = "character",
    help = "Name of marker by cell type reference file used with CellAssign.
            File name is expected to be in form: `<ref_name>_<source>_<version>.tsv`"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

# Set up -----------------------------------------------------------------------

#' Get reference information from a reference file
#'
#' @param ref_filename reference file name to parse
#' @extension file extension to consider during parsing
#'
#' @return list of reference file components to include in SCE metadata
get_ref_info <- function(ref_filename, extension) {
  ref_info <- ref_filename |>
    basename() |>
    # select everything before the extension
    stringr::word(1, sep = extension) |>
    # create a vector with name, source, version
    stringr::str_split(pattern = "_") |>
    unlist() |>
    purrr::set_names(c("ref_name", "ref_source", "ref_version"))

  return(ref_info)
}

# check that input file exists
if (!file.exists(opt$input_sce_file)) {
  stop("Missing input SCE file")
}

# check that output file ends in rds
if (!(stringr::str_ends(opt$output_sce_file, ".rds"))) {
  stop("output sce file name must end in .rds")
}

# read in input files
sce <- readr::read_rds(opt$input_sce_file)

# SingleR results --------------------------------------------------------------

if (!is.null(opt$singler_results)) {
  if (!file.exists(opt$singler_results)) {
    stop("Missing SingleR results file")
  }

  if (is.null(opt$singler_model_file)) {
    stop("SingleR model filename must be provided")
  }

  singler_results <- readr::read_rds(opt$singler_results)

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
  singler_barcodes <- annotations_df$barcodes
  sce_barcodes <- colnames(sce)
  missing_barcodes <- colnames(sce)[!(sce_barcodes %in% singler_barcodes)]

  # create a data frame with unclassified barcodes
  unclassified_df <- data.frame(
    barcodes = missing_barcodes,
    singler_celltype_annotation = "Unclassified cell"
  )

  # combine into one data frame with classified and unclassified cells
  all_annotations_df <- dplyr::bind_rows(list(annotations_df, unclassified_df))

  # add annotations to colData
  new_coldata <- colData(sce) |>
    as.data.frame() |>
    dplyr::left_join(all_annotations_df, by = c("barcodes")) |>
    DataFrame(row.names = colData(sce)$barcodes)
  colData(sce) <- new_coldata

  # get reference name, source and version
  singler_ref_info <- get_ref_info(
    ref_filename = opt$singler_model_file,
    extension = "_model\\.rds"
  )

  # add singler info to metadata
  metadata(sce)$singler_results <- singler_results
  metadata(sce)$singler_reference <- singler_ref_info[["ref_name"]]
  metadata(sce)$singler_reference_label <- label_type
  metadata(sce)$singler_reference_source <- singler_ref_info[["ref_source"]]
  metadata(sce)$singler_reference_version <- singler_ref_info[["ref_version"]]

  # add note about cell type method to metadata
  metadata(sce)$celltype_methods <- c(metadata(sce)$celltype_methods, "singler")
}

# CellAssign results -----------------------------------------------------------

if (!is.null(opt$cellassign_predictions)) {
  # check that cellassign predictions file was provided
  if (!file.exists(opt$cellassign_predictions)) {
    stop("Missing CellAssign predictions file")
  }

  if (is.null(opt$cellassign_ref_file)) {
    stop("CellAssign reference filename must be provided")
  }

  predictions <- readr::read_tsv(opt$cellassign_predictions)

  # get cell type with maximum prediction value for each cell
  celltype_assignments <- predictions |>
    tidyr::pivot_longer(
      !barcode,
      names_to = "celltype",
      values_to = "prediction"
    ) |>
    dplyr::group_by(barcode) |>
    dplyr::slice_max(prediction, n = 1) |>
    dplyr::ungroup()

  # join by barcode to make sure assignments are in the right order
  celltype_assignments <- data.frame(barcode = sce$barcodes) |>
    dplyr::left_join(celltype_assignments, by = "barcode") |>
    # any cells that are NA were not classified by cellassign
    dplyr::mutate(celltype = ifelse(is.na(celltype), "Unclassified cell", celltype))

  # add cell type and prediction to colData
  sce$cellassign_celltype_annotation <- celltype_assignments$celltype
  sce$cellassign_max_prediction <- celltype_assignments$prediction

  # get reference name, source and version
  cellassign_ref_info <- get_ref_info(
    ref_filename = opt$cellassign_ref_file,
    extension = "\\.tsv"
  )

  # add entire predictions matrix, ref name, and version to metadata
  metadata(sce)$cellassign_predictions <- predictions
  metadata(sce)$cellassign_reference <- cellassign_ref_info[["ref_name"]]
  metadata(sce)$cellassign_reference_source <- cellassign_ref_info[["ref_source"]]
  metadata(sce)$cellassign_reference_version <- cellassign_ref_info[["ref_version"]]

  # add cellassign as celltype method
  # note that if `metadata(sce)$celltype_methods` doesn't exist yet, this will
  #  come out to just the string "cellassign"
  metadata(sce)$celltype_methods <- c(metadata(sce)$celltype_methods, "cellassign")
}

# export annotated object with cellassign assignments
readr::write_rds(sce, opt$output_sce_file, compress = "gz")
