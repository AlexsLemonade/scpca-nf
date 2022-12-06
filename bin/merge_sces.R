#!/usr/bin/env Rscript

# import libraries
library(optparse)
library(SingleCellExperiment)

# Set up optparse options
option_list <- list(
  make_option(
    opt_str = c("--input_library_ids"),
    type = "character",
    help = "Comma separated list of library IDs corresponding to the libraries being integrated."
  ),
  make_option(
    opt_str = c("--input_sce_files"),
    type = "character",
    help = "Comma separated list of input sce file paths corresponding to the sces being integrated."
  ),
  make_option(
    opt_str = c("-o", "--output_sce_file"),
    type = "character",
    help = "Path to output RDS file, must end in .rds"
  )
)

# Setup ------------------------------------------------------------------------

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

# check that file extension for output file is correct 

# list of paths to sce files 
input_sce_files <- unlist(stringr::str_split(opt$input_sce_files, ','))

# set up library ids
if(is.null(opt$input_library_ids)){
  # extract library ids from filename
  input_library_ids <- stringr::word(basename(input_sce_files), 1, sep = "_")
} else {
  # pull library ids from list 
  input_library_ids <- unlist(stringr::str_split(opt$input_library_ids, ',')) 
}

# check that input files exist 
missing_sce_files <- input_sce_files[which(!file.exists(input_sce_files))]
if(length(missing_sce_files) > 0){
  stop(
    glue::glue(
      "\nMissing SCE object for {missing_sce_files}."
    )
  )
}
# Merge SCEs -------------------------------------------------------------------

# get list of sces
sce_list <- purrr::map(input_sce_files, readr::read_rds)
names(sce_list) <- input_library_ids

# create combined SCE object
merged_sce <- scpcaTools::merge_sce_list(sce_list,
                                         batch_column = "library_id",
                                         preserve_rowdata_cols = "gene_symbol",
                                         cell_id_column = "cell_id")


# HVG selection ----------------------------------------------------------------

# extract the column with the block variable
batch_column <- merged_sce$library_id

# model gene variance
gene_var_block <- scran::modelGeneVar(merged_sce,
                                      block = batch_column)
# identify subset of variable genes
hvg_list <- scran::getTopHVGs(gene_var_block,
                               n = 2000)

metadata(merged_sce)$merged_highly_variable_genes <- hvg_list

# Dim Reduction PCA and UMAP----------------------------------------------------

# multi batch PCA on merged object
multi_pca <- batchelor::multiBatchPCA(merged_sce,
                                      subset.row = hvg_list,
                                      batch = batch_column,
                                      preserve.single = TRUE)

# add PCA results to merged SCE object 
reducedDim(merged_sce, "PCA") <- multi_pca@listData[[1]]

# add UMAP 
merged_sce <- scater::runUMAP(dimred = "PCA")

# write out merged sce file 
readr::write_rds(merged_sce, opt$output_sce_file)
