#!/usr/bin/env Rscript

# This script is used to classify and annotate cells using SingleR

# import libraries
suppressPackageStartupMessages({
  library(optparse)
  library(SingleCellExperiment)
})

# set up arguments
option_list <- list(
  make_option(
    opt_str = c("-i", "--sce_file"),
    type = "character",
    help = "path to rds file with sce object to perform cell typing on"
  ),
  make_option(
    opt_str = c("--singler_model_file"),
    type = "character",
    help = "path to file containing a single model generated for SingleR annotation.
            File name is expected to be in form: <ref_name>_<source>_<version>_<gene_set_version>_model.rds."
  ),
  make_option(
    opt_str = c("--output_singler_results_file"),
    type = "character",
    help = "path to output RDS file that will store the SingleR results object. Must end in .rds"
  ),
  make_option(
    opt_str = c("--seed"),
    type = "integer",
    help = "A random seed for reproducibility."
  ),
  make_option(
    opt_str = c("-t", "--threads"),
    type = "integer",
    default = 1,
    help = "Number of multiprocessing threads to use."
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

# Set up -----------------------------------------------------------------------

# set seed
set.seed(opt$seed)

# check that input file file exists
if (!file.exists(opt$sce_file)) {
  stop("Missing SCE file")
}

# check that output files have the right extensions
if (!(stringr::str_ends(opt$output_singler_results_file, ".rds"))) {
  stop("output SingleR result file name must end in .rds")
}

# check that reference exists and filename is properly formatted
singler_model_file <- opt$singler_model_file
if (!file.exists(singler_model_file)) {
  stop(glue::glue("Provided model file {singler_model_file} is missing."))
}
if (!(stringr::str_ends(singler_model_file, "_model.rds"))) {
  stop(glue::glue("Provided model file {singler_model_file} must end in '_model.rds.'"))
}

# get & check reference name
reference_name <- stringr::str_remove(singler_model_file, "_model.rds$")
if (reference_name == "") {
  stop(glue::glue("Provided model file name must be formatted as `<model_name>_model.rds`"))
}

# set up multiprocessing params
if (opt$threads > 1) {
  bp_param <- BiocParallel::MulticoreParam(opt$threads)
} else {
  bp_param <- BiocParallel::SerialParam()
}

# read in input rds file
sce <- readr::read_rds(opt$sce_file)

# read in model
singler_model <- readr::read_rds(singler_model_file)

# perform cell typing annotation
singler_results <- SingleR::classifySingleR(
  trained = singler_model,
  test = sce,
  fine.tune = TRUE,
  BPPARAM = bp_param
)

# add reference name to singler_results DataFrame metadata
metadata(singler_results)$reference_name <- reference_name
# add label name to metadata
metadata(singler_results)$reference_label <- singler_model$reference_label
# save cell ontology table to results
# if this doesn't exist it will just be NULL
metadata(singler_results)$cell_ontology_df <- singler_model$cell_ontology_df

# export results ---------------

# next, the full result to a compressed rds
readr::write_rds(
  singler_results,
  opt$output_singler_results_file,
  compress = "bz2"
)
