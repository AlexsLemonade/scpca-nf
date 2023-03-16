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
    opt_str = c("--singler_models"),
    type = "character",
    help = "list of models generated for use with SingleR. Each input file contains 
      a list of models generated from a single reference, one each for each label type:
      `label.main`, `label.fine`, and `label.ont`."
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
model_files <- unlist(stringr::str_split(opt$singler_models, ","))
if(!all(file.exists(model_files))){
  missing_files <- model_files[which(!file.exists(model_files))]
  glue::glue("
             Missing model file(s): {missing_files}
             ")
  stop("Please make sure that all provided SingleR models exist.")
}

# set up multiprocessing params
if(opt$threads > 1){
  bp_param = BiocParallel::MulticoreParam(opt$threads)
} else {
  bp_param = BiocParallel::SerialParam()
}

# read in input rds file
sce <- readr::read_rds(opt$input_sce_file)

# read in references as a list of lists
# each file contains a named list of models generated using the same reference dataset
# but unique labels in the reference dataset
model_names <- stringr::str_remove(basename(model_files), "_model.rds")
names(model_files) <- model_names
model_list <- purrr::map(model_files, readr::read_rds) |>
  # ensure we have label type before reference name
  # example: label.main-HumanPrimaryCellAtlasData
  # where `label.main` is the name of the model stored in the file and
  # `HumanPrimaryCellAtlasData` is the name of the reference used for each file containing a list of models
  purrr::imap(\(model_list, ref_name){
                names(model_list) <- glue::glue("{names(model_list)}-{ref_name}")
                model_list
              }) |>
  purrr::flatten() 

# SingleR classify -------------------------------------------------------------

# create a partial function for mapping easily
classify_sce <- purrr::partial(SingleR::classifySingleR, 
                               test = sce, 
                               fine.tune=TRUE, 
                               BPPARAM = bp_param)
# run singleR for all provided models
all_singler_results <- model_list |>
    purrr::map(classify_sce)

# Annotate sce -----------------------------------------------------------------

# create a dataframe with a single column of annotations for each model used
all_annotations_df <- all_singler_results |>
  purrr::map_dfc(\(result) result$pruned.labels ) |>
  DataFrame(check.names = FALSE) # prevent replacing "-" in annotation columns

colData(sce) <- cbind(colData(sce), all_annotations_df)

# store results in metadata
metadata(sce)$singler_results <- all_singler_results


# export sce with annotations added
readr::write_rds(sce,
                 opt$output_sce_file,
                 compress = 'gz')

