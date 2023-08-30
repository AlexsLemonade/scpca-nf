#!/usr/bin/env Rscript

# This script is used to read in the predictions from CellAssign and assign cell types in the annotated RDS file

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
    opt_str = c("--cellassign_predictions"),
    type = "character",
    help = "path to tsv file containing the prediction matrix returned by running CellAssign"
  ),
  make_option(
    opt_str = c("--reference_name"),
    type = "character",
    help = "name referring to the marker gene reference used for CellAssign"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

# check that input file exists
if (!file.exists(opt$input_sce_file)){
  stop("Missing input SCE file")
}

# check that cellassign predictions file was provided
if (!file.exists(opt$cellassign_predictions)){
  stop("Missing CellAssign predictions file")
}
# check that reference_name was provided
if (is.null(opt$reference_name)) {
  stop("Missing reference name")
}

# check that output file ends in rds
if(!(stringr::str_ends(opt$output_sce_file, ".rds"))){
  stop("output sce file name must end in .rds")
}
# read in input files
sce <- readr::read_rds(opt$input_sce_file)
predictions <- readr::read_tsv(opt$cellassign_predictions)

celltype_assignments <- predictions |>
  tidyr::pivot_longer(!barcode,
                      names_to = "celltype",
                      values_to = "prediction") |>
  dplyr::group_by(barcode) |>
  dplyr::slice_max(prediction, n = 1) |>
  dplyr::ungroup()

# join by barcode to make sure assignments are in the right order
celltype_assignments <- data.frame(barcode = sce$barcodes) |>
  dplyr::left_join(celltype_assignments, by = "barcode")

# add cell type and prediction to colData
sce$cellassign_celltype_annotation <- celltype_assignments$celltype
sce$cellassign_max_prediction <- celltype_assignments$prediction

# add entire predictions matrix and ref name to metadata
metadata(sce)$cellassign_predictions <- predictions
metadata(sce)$cellassign_reference <- opt$reference_name

# export annotated object with cellassign assignments 
readr::write_rds(sce, opt$output_sce_file)
