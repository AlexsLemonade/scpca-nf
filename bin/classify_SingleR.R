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
    opt_str = c("-i", "--input_sce_file"),
    type = "character",
    help = "path to rds file with input sce object"
  ),
  make_option(
    opt_str = c("-o", "--output_sce_file"),
    type = "character",
    help = "path to output rds file to store processed sce object. Must end in .rds"
  ),
  make_option(
    opt_str = c("--singler_model_file"),
    type = "character",
    help = "path to file containing a single model generated for SingleR annotation"
  ),
  make_option(
    opt_str = c("--label_name"),
    type = "character",
    default = "label.ont",
    help = "label used when building the SingleR reference"
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
if(!file.exists(opt$input_sce_file)){
  stop("Missing input SCE file")
}

# check that references all exist
singler_model_file <- opt$singler_model_file
if(!file.exists(singler_model_file)) {
  stop(glue::glue("Provided model file {singler_model_file} is missing."))
}

# set up multiprocessing params
if(opt$threads > 1){
  bp_param = BiocParallel::MulticoreParam(opt$threads)
} else {
  bp_param = BiocParallel::SerialParam()
}

# read in input rds file
sce <- readr::read_rds(opt$input_sce_file)

# read in model
singler_model <- readr::read_rds(singler_model_file)

# SingleR classify and annotate SCE-------------------------------------------------------------

singler_results <- SingleR::classifySingleR(
  trained = singler_model,
  test = sce,
  fine.tune = TRUE,
  BPPARAM = bp_param
)

# add annotations to SCE colData
if (opt$label_name == "label.ont") {
  # If ontologies were used, create ontology column and join in the cell names
  sce$singler_celltype_ontology <- singler_results$pruned.labels

  colData(sce) <- colData(sce) |>
    as.data.frame() |>
    dplyr::left_join(
      singler_model$cell_ontology_df,
      by = c("singler_celltype_ontology" = "ontology_id")
    ) |>
    # this is the new column that was joined in with the cell names
    dplyr::rename(singler_celltype_annotation = ontology_cell_names) |>
    # make sure we keep rownames
    DataFrame(row.names = colData(sce)$barcodes)


} else {
  # otherwise, just add cell names
  sce$singler_celltype_annotation <- singler_results$pruned.labels
}

# store full SingleR results in metadata, as well as other SingleR parameters
metadata(sce)$singler_results <- singler_results
metadata(sce)$singler_reference <- singler_model$reference_name
metadata(sce)$singler_reference_label <- singler_model$reference_label


# export sce with annotations added
readr::write_rds(sce,
                 opt$output_sce_file,
                 compress = 'gz')

