#!/usr/bin/env Rscript

# This script is used to build and train a model to be used for SingleR.
# The expected input is a SummarizedExperiment containing a reference dataset and a
# 3-column, transcript to gene tsv file (as used with Alevin-fry), to provide the expected
# geneset to be used for training the reference dataset. The 3 columns in this file should
# contain the transcript id, gene id, and the transcript type in that order, but need not
# contain column names.

# import libraries
suppressPackageStartupMessages({
  library(optparse)
  library(SingleCellExperiment)
})

# set up arguments
option_list <- list(
  make_option(
    opt_str = c("--ref_file"),
    type = "character",
    help = "path to rds file with reference dataset to use for cell type annotation.
      These reference datasets must contain annotations labeld with 
      `label.main`, `label.fine`, and `label.ont`."
  ),
  make_option(
    opt_str = c("--output_file"),
    type = "character",
    help = "path to output rds file to store trained model for reference dataset."
  ),
  make_option(
    opt_str = c("--fry_tx2gene"),
    type = "character",
    help = "path to tsv file containing three columns with transcript id, gene id, and
      type of transcript (either spliced or unspliced)."
  ),
  make_option(
    opt_str = c("-r", "--random_seed"),
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

# check that input files exist
if(!file.exists(opt$ref_file)){
  stop("Missing input file with cell type reference (`ref_file`).")
}

if(!file.exists(opt$fry_tx2gene)){
  stop("Missing `fry_tx2gene` file.")
}

# set up multiprocessing params
if(opt$threads > 1){
  bp_param = BiocParallel::MulticoreParam(opt$threads)
} else {
  bp_param = BiocParallel::SerialParam()
}

# read in model
ref_data <- readr::read_rds(opt$ref_file)

# check that ref data contains correct labels
if(!all(c("label.main", "label.fine", "label.ont") %in% colnames(colData(ref_data)))){
  stop("Reference dataset must contain `label.main`, `label.fine`, and `label.ont` in `colData`.")
}

# read in tx2gene
tx2gene <- readr::read_tsv(opt$fry_tx2gene,
                           col_names = c("transcript", "gene", "transcript_type"))

# select genes to use for model restriction
gene_ids <- unique(tx2gene$gene)

# check that genes aren't empty
if (length(gene_ids) == 0){
  stop("Provided tx2gene tsv file does not contain any genes.")
}

# Train models -----------------------------------------------------------------

label_cols <- c("label.main", "label.fine", "label.ont")
names(label_cols) <- label_cols

models <- label_cols |>
 purrr::map(\(label_col) {
   SingleR::trainSingleR(
     ref_data,
     labels = colData(ref_data)[, label],
     genes = "de",
     # only use genes found in index
     restrict = gene_ids,
     BPPARAM = bp_param    
   )})

# export models
readr::write_rds(models, opt$output_file)