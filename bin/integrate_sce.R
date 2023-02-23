#!/usr/bin/env Rscript

# Script used to perform integration on a given merged SCE object using R-based methods
#
# A merged SCE file in RDS format is read in. Integration is performed with either
# `fastMNN` or `harmony`. The integrated SCE object is saved as an RDS
# file in the provided output file.
#

# import libraries
library(optparse)
library(SingleCellExperiment)

# Set up optparse options
option_list <- list(
  make_option(
    opt_str = c("-i", "--input_sce_file"),
    type = "character",
    default = NULL,
    help = "Path to RDS file that contains the merged SCE object to integrate"
  ),
  make_option(
    opt_str = c("-o", "--output_sce_file"),
    type = "character",
    default = NULL,
    help = "Path to RDS file where the integrated SCE object will be saved"
  ),
  make_option(
    opt_str = c("--method"),
    type = "character",
    default = NULL,
    help = "Integration method to use, either `fastMNN` or `harmony` (case-sensitive)."
  ),
  make_option(
    opt_str = c("--harmony_covariate_cols"),
    type = "character",
    default = NULL,
    help = "Optional comma-separated list of columns (e.g. patient, sex) to consider as covariates
            during integration with `harmony`."
  ),
  make_option(
    opt_str = c("-t", "--threads"),
    type = "integer",
    default = 1,
    help = "Number of multiprocessing threads to use"
  ),
  make_option(
    opt_str = c("--seed"),
    type = "integer",
    default = NULL,
    help = "random seed to set during integration"
  )
)

# Setup ------------------------------------------------------------------------
# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

set.seed(opt$seed)

# Check and assign provided method based on available methods
available_methods <- c("fastMNN", "harmony")

if(!(opt$method %in% available_methods)){
  stop("You must specify either `fastMNN` or `harmony` to the --method.")
}

# Check that provided input file exists and is an RDS file
if(is.null(opt$input_sce_file)) {
  stop("You must provide the path to the RDS file with merged SCEs to --input_sce_file")
} else {
  if(!file.exists(opt$input_sce_file)) {
    stop("Provided --input_sce_file file does not exist.")
  }
}

# set up multiprocessing params
if(opt$threads > 1){
  bp_param = BiocParallel::MulticoreParam(opt$threads)
} else {
  bp_param = BiocParallel::SerialParam()
}

# Read in SCE file -------------------------------------------------------------
merged_sce <- readr::read_rds(opt$input_sce_file)

# check that input contains a SCE object
if(!is(merged_sce, "SingleCellExperiment")){
  stop("The input RDS file must contain a SingleCellExperiment object.")
}

# hvgs to use for integration
merged_hvgs <- metadata(merged_sce)$merged_hvgs

# Integration ------------------------------------------------------------------

# create list of arguments for either integration method
integration_args <- list(
  merged_sce = merged_sce,
  integration_method = opt$method,
  batch_column = "library_id",
  seed = opt$seed,
  return_corrected_expression = FALSE
)

# append fastMNN only args
if(opt$method == "fastMNN"){
  integration_args <-  append(
    integration_args,
    list(
      subset.row = merged_hvgs,
      auto.merge = TRUE,
      BPPARAM = bp_param
   ))
}

# append harmony only args
if(opt$method == "harmony"){
  integration_arts <- append(
    integration_args,
    list(
      covariate_cols = opt$harmony_covariate_cols
    ))
}

# perform integration
integrated_sce <- do.call(scpcaTools::integrate_sces, integration_args)

# calculate UMAP from corrected PCA
integrated_sce <- integrated_sce |>
  scater::runUMAP(dimred = glue::glue("{opt$method}_PCA"),
                  name = glue::glue("{opt$method}_UMAP"))

# add method to integrated sce
existing_methods <- metadata(integrated_sce)$integration_methods
if(!is.null(existing_methods)){
  metadata(integrated_sce)$integration_methods <- c(existing_methods, opt$method)
} else {
  metadata(integrated_sce)$integration_methods <- opt$method
}

# write out integrated object with merged data + corrected PCA and UMAP
readr::write_rds(integrated_sce, opt$output_sce_file)
