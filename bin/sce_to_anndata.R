#!/usr/bin/env Rscript

# This script takes a SingleCellExperiment stored in a .rds file and converts the main experiment
# (usually RNA) to an AnnData object saved as an hdf5 file

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
    help = "path to rds file with input sce object to be converted"
  ),
  make_option(
    opt_str = c("--output_rna_h5"),
    type = "character",
    help = "path to output hdf5 file to store RNA counts as AnnData object. Must end in .hdf5 or .h5"
  ),
  make_option(
    opt_str = c("--output_feat_h5"),
    type = "character",
    help = "path to output hdf5 file to store feature counts as AnnData object. 
    Only used if the input SCE contains an altExp. Must end in .hdf5 or .h5"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

# Set up -----------------------------------------------------------------------

# check that filtered SCE file exists
if(!file.exists(opt$input_sce_file)){
  stop(glue::glue("{opt$input_sce_file} does not exist."))
}

# check that output file is h5
if(!(stringr::str_ends(opt$output_rna_h5, ".hdf5|.h5"))){
  stop("output rna file name must end in .hdf5 or .h5")
}

# Convert to AnnData -----------------------------------------------------------

# read in sce
sce <- readr::read_rds(opt$input_sce_file)

# export sce as anndata object
scpcaTools::sce_to_anndata(
  sce,
  anndata_file = opt$output_rna_h5
)

# check for altExp
alt_name <- altExpNames(sce)

# if altExp exists then export as hdf5
if(!is.null(alt_name) & length(alt_name) == 1){
  
  # check for output file 
  if(!(stringr::str_ends(opt$output_feat_h5, ".hdf5|.h5"))){
    stop("output feature file name must end in .hdf5 or .h5")
  }
  
  # extract altExp 
  alt_sce <- altExp(sce, alt_name)
  
  # export altExp sce as anndata object
  scpcaTools::sce_to_anndata(
    alt_sce,
    anndata_file = opt$output_feat_h5
  )
  
}

# warn that only one altExp will be exported 
if (!is.null(alt_name) & length(alt_name) > 1){
  warning("Only 1 altExp named with {alt_name[1]} will be exported to HDF5")
}
