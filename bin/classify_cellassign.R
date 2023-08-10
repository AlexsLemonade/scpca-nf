#!/usr/bin/env Rscript

# This script is used to read in the predictions from CellAssign and assign cell types

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

# check that input file file exists
if(!file.exists(opt$input_sce_file)){
  stop("Missing input SCE file")
}

if(!file.exists(opt$cellassign_predictions)){
  stop("Missing CellAssign predictions file")
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

sce$cellassign_celltype_annotation <- celltype_assignments

metadata(sce)$cellassign_predictions <- predictions
metadata(sce)$cellassign_reference <- opt$reference_name
