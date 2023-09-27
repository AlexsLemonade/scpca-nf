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
    help = "path to file containing a single model generated for SingleR annotation"
  ),
  make_option(
    opt_str = c("--output_singler_annotations_file"),
    type = "character",
    help = "path to output TSV file that will store the SingleR annotations. Must end in .rds"
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
set.seed(opt$random_seed)

# check that input file file exists
if (!file.exists(opt$sce_file)) {
  stop("Missing SCE file")
}

# check that output files have the righr extensions
if (!(stringr::str_ends(opt$output_singler_results_file, ".rds"))) {
  stop("output SingleR result file name must end in .rds")
}
if (!(stringr::str_ends(opt$output_singler_annotations_file, ".tsv"))) {
  stop("output SingleR annotations file name must end in .tsv")
}

# check that references all exist
singler_model_file <- opt$singler_model_file
if (!file.exists(singler_model_file)) {
  stop(glue::glue("Provided model file {singler_model_file} is missing."))
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

# export results

# first, a stand-alone tsv of annotations with both pruned and full labels
readr::write_tsv(
  tibble::tibble(
    barcode = rownames(singler_results),
    pruned_labels = singler_results$pruned_labels
  ),
  output_singler_annotations_file
)

# next, the full result to a compressed rds
readr::write_rds(
  singler_results,
  opt$output_singler_results_file,
  compress = "gz"
)
