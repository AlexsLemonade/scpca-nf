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
    opt_str = c("--feature_name"),
    type = "character",
    help = "Feature type. Must match the altExp name, if present."
  ),
  make_option(
    opt_str = c("--output_feature_h5"),
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

# grab sample metadata
sample_metadata <- metadata(sce)$sample_metadata

# add library id as a column to the sce object
sce$library_id <- metadata(sce)$library_id

# add sample metadata to colData sce
sce <- scpcaTools::metadata_to_coldata(sce,
                                       join_columns = "library_id")
# remove sample metadata from sce metadata, otherwise conflicts with converting object
metadata(sce) <- metadata(sce)[names(metadata(sce)) != "sample_metadata"]

# export sce as anndata object
scpcaTools::sce_to_anndata(
  sce,
  anndata_file = opt$output_rna_h5
)

# if feature data exists, grab it and export to AnnData
if(!is.null(opt$feature_name)){

  # make sure the feature data is present
  if(!(opt$feature_name %in% altExpNames(sce))){
    stop("feature_name must match name of altExp in provided SCE object.")
  }

  # check for output file
  if(!(stringr::str_ends(opt$output_feature_h5, ".hdf5|.h5"))){
    stop("output feature file name must end in .hdf5 or .h5")
  }

  # extract altExp
  alt_sce <- altExp(sce, opt$feature_name)

  # add library ID to colData
  alt_sce$library_id <- metadata(sce)$library_id

  # add sample metadata to alt sce metadata
  metadata(alt_sce)$sample_metadata <- sample_metadata

  # add sample metadata to alt sce coldata
  alt_sce <- scpcaTools::metadata_to_coldata(alt_sce,
                                             join_columns = "library_id")

  # remove sample metadata from metadata, otherwise conflicts with converting object
  metadata(sce) <- metadata(sce)[!metadata(sce) %in% sample_metadata]

  # export altExp sce as anndata object
  scpcaTools::sce_to_anndata(
    alt_sce,
    anndata_file = opt$output_feature_h5
  )

}
