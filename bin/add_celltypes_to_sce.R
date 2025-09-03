#!/usr/bin/env Rscript

# This script is used to read in and assign these cell types in the annotated RDS file:
# - CellAssign
# - SingleR
# - consensus cell types
# This script additionally counts the number of normal reference cells in the library

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
    help = "path to rds file with input sce object",
    default = ""
  ),
  make_option(
    opt_str = c("-o", "--output_sce_file"),
    type = "character",
    help = "path to output rds file to store annotated sce object. Must end in .rds"
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
    opt_str = c("--consensus_celltype_ref"),
    type = "character",
    help = "Path to file containing the reference for assigning consensus cell type labels",
    default = ""
  ),
  make_option(
    opt_str = c("--consensus_validation_ref"),
    type = "character",
    help = "Path to file mapping consensus validation groups to consensus labels for counting normal reference cells intended as input to inferCNV",
    default = ""
  ),
  make_option(
    opt_str = c("--diagnosis_celltype_ref"),
    type = "character",
    help = "Path to file mapping broad diagnoses to consensus validation groups for counting normal reference cells intended as input to inferCNV",
    default = ""
  ),
  make_option(
    opt_str = c("--diagnosis_groups_ref"),
    type = "character",
    help = "Path to file mapping broad diagnoses to individual diagnoses for counting normal reference cells intended as input to inferCNV",
    default = ""
  ),
  make_option(
    opt_str = c("--normal_cells_file"),
    type = "character",
    help = "Path to write number of calculated normal cells to.
      This calculation is only performed if `diagnosis_celltype_ref` and `diagnosis_groups_ref` are provided.
      If not calculated, this file will be created as an empty file",
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
    names(ref_info) <- c("ref_name", "ref_source", "ref_version", "gene_set_version", "date")
  } else if (ref_type == "CellAssign") {
    names(ref_info) <- c("ref_name", "ref_source", "ref_version")
  }

  return(ref_info)
}

# check that input file exists and output file ends in rds
stopifnot(
  "Missing input SCE file" = file.exists(opt$input_sce_file),
  "output sce file name must end in .rds" = stringr::str_ends(opt$output_sce_file, ".rds"),
  "output file to store counted normal cells was not provided" = !is.null(opt$normal_cells_file)
)

# read in input files
sce <- readr::read_rds(opt$input_sce_file)

# We'll count normal cells when assigning consensus cell types
# Set to an empty string by default, to write out if we don't have consensus
normal_cells <- ""

# SingleR results --------------------------------------------------------------

has_singler <- file.exists(opt$singler_results)
if (has_singler) {
  # check singler model has been provided
  stopifnot("Singler model filename must be provided" = opt$singler_model_file != "")

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
    ref_filename = opt$singler_model_file,
    extension = "_model\\.rds",
    ref_type = "SingleR"
  )

  # add singler info to metadata
  metadata(sce)$singler_results <- singler_results
  metadata(sce)$singler_reference <- singler_ref_info[["ref_name"]]
  metadata(sce)$singler_reference_label <- label_type
  metadata(sce)$singler_reference_source <- singler_ref_info[["ref_source"]]
  metadata(sce)$singler_reference_version <- singler_ref_info[["ref_version"]]
  metadata(sce)$singler_gene_set_version <- singler_ref_info[["gene_set_version"]]
  metadata(sce)$singler_date <- singler_ref_info[["date"]]

  # add note about cell type method to metadata
  metadata(sce)$celltype_methods <- c(metadata(sce)$celltype_methods, "singler")
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

  # read in panglao ontology reference
  panglao_ref_df <- readr::read_tsv(opt$panglao_ontology_ref)

  # if cell assign predictions file exists but is empty then cell assign failed and we want to account for that with Not Run
  if (file.size(opt$cellassign_predictions) > 0) {
    # read in predictions file
    predictions <- readr::read_tsv(opt$cellassign_predictions)
  } else {
    predictions <- NULL
    has_cellassign <- FALSE # reset to false so that we don't add in consensus cell types
  }

  # if the only column is the barcode column or if the predictions file was empty
  # then CellAssign didn't complete successfully
  # otherwise add in cell type annotations and metadata to SCE
  if (is.null(predictions) || all(colnames(predictions) == "barcode")) {
    # if failed then note that in the cell type column
    sce$cellassign_celltype_annotation <- "Not run"
    has_cellassign <- FALSE # reset to false if cellassign didn't complete so we don't add consensus
  } else {
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
      dplyr::mutate(celltype = ifelse(is.na(celltype), "Unclassified cell", celltype))

    # add cell type and prediction to colData
    sce$cellassign_celltype_annotation <- celltype_assignments$celltype
    sce$cellassign_celltype_ontology <- celltype_assignments$ontology
    sce$cellassign_max_prediction <- celltype_assignments$prediction

    # get reference name, source and version
    cellassign_ref_info <- get_ref_info(
      ref_filename = opt$cellassign_ref_file,
      extension = "\\.tsv",
      ref_type = "CellAssign"
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

    # add cellassign reference organs to metadata
    cellassign_organs <- opt$celltype_ref_metafile |>
      readr::read_tsv() |>
      dplyr::filter(celltype_ref_name == cellassign_ref_info[["ref_name"]]) |>
      dplyr::pull(organs)

    if (cellassign_organs == "" | is.na(cellassign_organs)) {
      stop("Failed to obtain CellAssign reference organ list.")
    }
    metadata(sce)$cellassign_reference_organs <- cellassign_organs
  }
}

# assign consensus cell type labels
if (has_singler && has_cellassign) {
  # now make sure that reference file exists
  stopifnot(
    "Consensus cell type reference file does not exist" = file.exists(opt$consensus_celltype_ref)
  )

  # read in consensus table
  consensus_ref_df <- readr::read_tsv(opt$consensus_celltype_ref) |>
    # select unique combinations of consensus refs based on ontology columns
    dplyr::select(blueprint_ontology, panglao_ontology, consensus_ontology, consensus_annotation) |>
    dplyr::distinct()

  # create df with consensus assignments
  celltype_df <- colData(sce) |>
    as.data.frame() |>
    dplyr::select(
      barcodes,
      contains("celltype") # get both singler and cellassign with ontology
    ) |>
    # then add consensus labels
    dplyr::left_join(
      consensus_ref_df,
      by = c(
        "singler_celltype_ontology" = "blueprint_ontology",
        "cellassign_celltype_ontology" = "panglao_ontology"
      ),
      relationship = "many-to-many" # account for multiple of the same cell type
    ) |>
    # use unknown for NA annotation but keep ontology ID as NA
    tidyr::replace_na(list(consensus_annotation = "Unknown"))

  # add consensus cell type and ontology to sce
  sce$consensus_celltype_annotation <- celltype_df$consensus_annotation
  sce$consensus_celltype_ontology <- celltype_df$consensus_ontology


  # Count the number of normal reference cells -------------------------

  # Recalculate if these optional reference files were provided
  # they have no defaults, so they will be NULL if not provided
  check_input_files <- c(
    opt$consensus_validation_ref,
    opt$diagnosis_celltype_ref,
    opt$diagnosis_groups_ref
  )

  # Only proceed if _none_ are empty strings
  if (all(check_input_files != "")) {
    stopifnot(
      "Validation cell type reference file does not exist" = file.exists(opt$consensus_validation_ref),
      "Diagnosis/cell type map reference file does not exist" = file.exists(opt$diagnosis_celltype_ref),
      "Diagnosis map reference file does not exist" = file.exists(opt$diagnosis_groups_ref)
    )

    consensus_validation_df <- readr::read_tsv(opt$consensus_validation_ref)
    diagnosis_celltype_df <- readr::read_tsv(opt$diagnosis_celltype_ref)
    diagnosis_map_df <- readr::read_tsv(opt$diagnosis_groups_ref)

    sample_diagnosis <- metadata(sce)$sample_metadata$diagnosis

    diagnosis_groups <- diagnosis_map_df |>
      # get the broad diagnosis
      dplyr::filter(submitted_diagnosis == sample_diagnosis) |>
      dplyr::pull("diagnosis_group")

    reference_celltypes <- diagnosis_celltype_df |>
      # get the cell type groups to consider for this diagnosis
      dplyr::filter(diagnosis_group %in% diagnosis_groups)
    dplyr::select(celltype_groups) |>
      tidyr::separate_rows(celltype_groups, sep = ", ") |>
      # get the consensus cell types
      dplyr::left_join(consensus_validation_df, by = c("celltype_groups" = "validation_group_annotation")) |>
      dplyr::pull(consensus_annotation) |>
      unique()

    normal_cells <- sum(celltype_df$consensus_annotation %in% reference_celltypes)
  }
}

# export annotated object with cell type assignments
readr::write_rds(sce, opt$output_sce_file, compress = "bz2")

# export normal cell count
readr::write_lines(normal_cells, opt$normal_cells_file)
