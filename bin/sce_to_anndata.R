#!/usr/bin/env Rscript

# This script takes a SingleCellExperiment stored in a .rds file and converts the main experiment
# (usually RNA) to an AnnData object saved as an hdf5 file

# The AnnData object being exported by this script is formatted to fit CZI schema: 3.0.0

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
    Only used if the input SCE contains an altExp. Must end in .hdf5, .h5, or .h5ad"
  ),
  make_option(
    opt_str = c("--compress_output"),
    action = "store_true",
    default = FALSE,
    help = "Compress the HDF5 file containing the AnnData object"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

# Set up -----------------------------------------------------------------------

# check that filtered SCE file exists
if (!file.exists(opt$input_sce_file)) {
  stop(glue::glue("{opt$input_sce_file} does not exist."))
}

# check that output file is h5
if (!(stringr::str_ends(opt$output_rna_h5, ".hdf5|.h5|.h5ad"))) {
  stop("output rna file name must end in .hdf5, .h5, or .h5ad")
}

# CZI compliance function ------------------------------------------------------

# this function applies any necessary reformatting or changes needed to make
# sure that the sce that is getting converted to AnnData is compliant with
# CZI 3.0.0 requirements: https://github.com/chanzuckerberg/single-cell-curation/blob/b641130fe53b8163e50c39af09ee3fcaa14c5ea7/schema/3.0.0/schema.md
format_czi <- function(sce) {
  # add schema version
  metadata(sce)$schema_version <- "3.0.0"

  # add library_id as an sce colData column
  # need this column to join in the sample metadata with the colData
  if (!("library_id" %in% colnames(colData(sce)))) {
    sce$library_id <- metadata(sce)$library_id
  }

  # add sample metadata to colData sce
  sce <- scpcaTools::metadata_to_coldata(
    sce,
    join_columns = "library_id"
  )

  # modify colData to be AnnData and CZI compliant
  coldata_df <- colData(sce) |>
    as.data.frame() |>
    dplyr::mutate(
      # create columns for assay and suspension ontology terms
      assay_ontology_term_id = metadata(sce)$assay_ontology_term_id,
      suspension_type = metadata(sce)$seq_unit,
      # add is_primary_data column; only needed for anndata objects
      is_primary_data = FALSE
    )

  # add colData back to sce object
  colData(sce) <- DataFrame(
    coldata_df,
    row.names = rownames(colData(sce))
  )

  # remove sample metadata from sce metadata, otherwise conflicts with converting object
  metadata(sce) <- metadata(sce)[names(metadata(sce)) != "sample_metadata"]

  # modify rowData
  # we don't do any gene filtering between normalized and raw counts matrix
  # so everything gets set to false
  rowData(sce)$feature_is_filtered <- FALSE

  # paste X to reduced dim names if present
  if (!is.null(reducedDimNames(sce))) {
    reducedDimNames(sce) <- glue::glue("X_{reducedDimNames(sce)}")
  }

  return(sce)
}

# MainExp to AnnData -----------------------------------------------------------

# read in sce
sce <- readr::read_rds(opt$input_sce_file)

# grab sample metadata
# we need this if we have any feature data that we need to add it o
sample_metadata <- metadata(sce)$sample_metadata

# make main sce czi compliant
sce <- format_czi(sce)

# export sce as anndata object
# this function will also remove any R-specific object types from the SCE metadata
#   before converting to AnnData
scpcaTools::sce_to_anndata(
  sce,
  anndata_file = opt$output_rna_h5,
  compression = ifelse(opt$compress_output, "gzip", "none")
)

# AltExp to AnnData -----------------------------------------------------------

# if feature data exists, grab it and export to AnnData
if (!is.null(opt$feature_name)) {
  # make sure the feature data is present
  if (!(opt$feature_name %in% altExpNames(sce))) {
    stop("feature_name must match name of altExp in provided SCE object.")
  }

  # check for output file
  if (!(stringr::str_ends(opt$output_feature_h5, ".hdf5|.h5|.h5ad"))) {
    stop("output feature file name must end in .hdf5, .h5, or .h5ad")
  }

  # extract altExp
  alt_sce <- altExp(sce, opt$feature_name)

  # only convert altExp with > 1 rows
  if (nrow(alt_sce) > 1) {
    # add sample metadata from main sce to alt sce metadata
    metadata(alt_sce)$sample_metadata <- sample_metadata

    # make sce czi compliant
    alt_sce <- format_czi(alt_sce)

    # export altExp sce as anndata object
    scpcaTools::sce_to_anndata(
      alt_sce,
      anndata_file = opt$output_feature_h5,
      compression = ifelse(opt$compress_output, "gzip", "none")
    )
  } else {
    # warn that the altExp cannot be converted
    message(
      glue::glue("
        Only 1 row found in altExp named: {opt$feature_name}.
        This altExp will not be converted to an AnnData object.
      ")
    )
  }
}
