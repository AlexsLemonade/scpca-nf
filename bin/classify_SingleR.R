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
    help = "list of models generated for use with SingleR"
  ),
  make_option(
    opt_str = c("-seed"),
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
model_files <- opt$singler_models
if(!any(file.exists(model_files))){
  missing_files <- model_files[!which(file.exists(model_files))]
  glue::glue("
             Missing model file: {missing_files}
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

# read in references as a list
model_names <- stringr::str_remove(basename(model_files), "_model.rds")
names(model_files) <- model_names
model_list <- purrr::map(model_files, readr::read_rds) |>
  # make sure that the labels have label type before reference dataset name
  purrr::map2(model_names, 
              \(model, name){
    names(model) <- glue::glue("{names(model)}-{name}")
    model
  }) |> 
  unname() |> # remove existing package names
  unlist() # need to collapse into one list of all models

# SingleR classify -------------------------------------------------------------

# run singleR for all provided models 
all_singler_results <- purrr::map(model_list,
                                  \(model) SingleR::classifySingleR(
                                    sce,
                                    model,
                                    fine.tune = TRUE,
                                    BPPARAM = bp_param
                                  )) |> 
  purrr::set_names(names(model_list))

# Annotate sce -----------------------------------------------------------------

# create a dataframe with a single column of annotations for each model used 
all_annotations_df <- purrr::map(all_singler_results, 
                                 \(result){result$pruned.labels}) |> 
  dplyr::bind_cols() |>
  dplyr::mutate(barcode = colnames(sce))

# extract coldata from sce object to join with annotations
coldata_df <- as.data.frame(colData(sce)) |>
  tibble::rownames_to_column("barcode") |>
  dplyr::left_join(all_annotations_df)

# add coldata with annotations back to sce
colData(sce) <- DataFrame(coldata_df, 
                          row.names = coldata_df$barcode, 
                          # prevent replacing "-" in annotation columns
                          check.names = FALSE)

# add in scores and delta to metadata of sce object
metadata(sce)$singler_score <- purrr::map(all_singler_results,
                                          \(result){result$scores})

metadata(sce)$singler_delta <- purrr::map(all_singler_results,
                                          \(result) {result$delta.next})


# export sce with annotations added
readr::write_rds(sce, opt$output_sce_file)

