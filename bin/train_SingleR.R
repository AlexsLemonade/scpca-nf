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
      These reference datasets must contain annotations labeled with the given
      `label_name` (default 'label.ont')."
  ),
  make_option(
    opt_str = c("--output_file"),
    type = "character",
    help = "path to output rds file to store trained model for reference dataset."
  ),
  make_option(
    opt_str = c("--label_name"),
    type = "character",
    default = "label.ont",
    help = "name of cell type label to use for training the reference dataset."
  ),
  make_option(
    opt_str = c("--fry_tx2gene"),
    type = "character",
    default = "",
    help = "path to tsv file containing three columns with transcript id, gene id, and
      type of transcript (either spliced or unspliced)."
  ),
  make_option(
    opt_str = c("--flex_probeset"),
    type = "character",
    default = "",
    help = "path to csv file with gene ids for probes used with 10x flex platform"
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
stopifnot(
  "Missing input file with cell type reference (`ref_file`)." = file.exists(opt$ref_file),
  "Either fry_tx2gene or flex_probeset must be provided" = file.exists(opt$fry_tx2gene) | file.exists(opt$flex_probeset)
)

# set up multiprocessing params
if (opt$threads > 1) {
  bp_param <- BiocParallel::MulticoreParam(opt$threads)
} else {
  bp_param <- BiocParallel::SerialParam()
}

# read in model
ref_data <- readr::read_rds(opt$ref_file)

# check that ref data contains correct labels, and set up label column for later
label_col <- opt$label_name
if (!label_col %in% colnames(colData(ref_data))) {
  stop(
    glue::glue("Reference dataset must contain `{label_col}` in `colData`.")
  )
}

if (file.exists(opt$fry_tx2gene)) {
  # read in tx2gene
  tx2gene <- readr::read_tsv(
    opt$fry_tx2gene,
    col_names = c("transcript", "gene", "transcript_type")
  )

  # select genes to use for model restriction
  gene_ids <- unique(tx2gene$gene)
} else if (file.exists(opt$flex_probeset)) {
  # read in probeset from 10x
  probes_df <- readr::read_csv(
    opt$flex_probeset,
    skip = 5
  )
  gene_ids <- unique(probes_df$gene_id)
}


# check that genes aren't empty
if (length(gene_ids) == 0) {
  stop("Provided tx2gene tsv file does not contain any genes.")
}

# Train model -----------------------------------------------------------------

singler_model <- SingleR::trainSingleR(
  ref_data,
  labels = colData(ref_data)[, label_col],
  genes = "de",
  # only use genes found in index
  restrict = gene_ids,
  BPPARAM = bp_param
)

# If the label is `label.ont` (default), save corresponding cell label names
if (label_col == "label.ont") {
  cl_ont <- ontoProc::getOnto("cellOnto")

  # grab ontology ids from singler_model
  ontology_ids <- singler_model$labels$unique

  # grab corresponding cell names
  cell_names <- cl_ont$name[ontology_ids]

  # save data frame with _matching_ ontology ids and cell names to model object
  # use tibble to avoid rownames
  singler_model$cell_ontology_df <- tibble::tibble(
    ontology_id = names(cell_names),
    ontology_cell_names = unname(cell_names)
  )
}

# Save reference name and label to model object
singler_model$reference_name <- stringr::str_replace(opt$output_file, "_model.rds$", "")
singler_model$reference_label <- label_col

# export models
readr::write_rds(singler_model, opt$output_file)
