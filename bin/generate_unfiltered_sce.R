#!/usr/bin/env Rscript

# This script takes the output folder from alevin-fry as input and
# returns the unfiltered counts matrices as a SingleCellExperiment stored in a .rds file

# import libraries
suppressPackageStartupMessages({
  library(optparse)
  library(SingleCellExperiment)
  library(scpcaTools)
})
# set up arguments
option_list <- list(
  make_option(
    opt_str = c("-a", "--alevin_dir"),
    type = "character",
    help = "directory with alevin output files for RNA-seq quantification"
  ),
  make_option(
    opt_str = c("-f", "--feature_dir"),
    type = "character",
    default = "",
    help = "directory with alevin output files for feature quantification"
  ),
  make_option(
    opt_str = c("-n", "--feature_name"),
    type = "character",
    default = "ALT",
    help = "Feature type"
  ),
  make_option(
    opt_str = c("-u", "--unfiltered_file"),
    type = "character",
    help = "path to output unfiltered rds file. Must end in .rds"
  ),
  make_option(
    opt_str = c("-m", "--mito_file"),
    type = "character",
    help = "path to list of mitochondrial genes"
  ),
  make_option(
    opt_str = c("-g", "--gtf_file"),
    type = "character",
    help = "path to gtf file with gene annotations"
  ),
  make_option(
    opt_str = c("-t", "--technology"),
    type = "character",
    help = "sequencing technology string to store in metadata"
  ),
  make_option(
    opt_str = c("--assay_ontology_term_id"),
    type = "character",
    default = NULL,
    help = "Experimental Factor Ontology term associated with provided tech_version"
  ),
  make_option(
    opt_str = c("-s", "--seq_unit"),
    type = "character",
    help = "sequencing unit string to store in metadata (e.g., cell, nucleus)"
  ),
  make_option(
    opt_str = c("--library_id"),
    type = "character",
    help = "library id"
  ),
  make_option(
    opt_str = c("--sample_id"),
    type = "character",
    help = "sample id(s). If more than one, separated by commas and/or semicolons."
  ),
  make_option(
    opt_str = c("--project_id"),
    type = "character",
    help = "project id"
  ),
  make_option(
    opt_str = c("--sample_metadata_file"),
    type = "character",
    help = "path to tsv file containing sample metadata"
  ),
  make_option(
    opt_str = c("--spliced_only"),
    action = "store_true",
    default = FALSE,
    help = "include only the spliced counts as the main counts assay in the returned SCE object"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

# check that output file name ends in .rds
if (!(stringr::str_ends(opt$unfiltered_file, ".rds"))) {
  stop("unfiltered file name must end in .rds")
}

# check that mitochondrial gene list exists
if (!file.exists(opt$mito_file)) {
  stop("Mitochondrial gene list file not found.")
}

# check that gtf file exists
if (!file.exists(opt$gtf_file)) {
  stop("gtf file not found.")
}

# check that sample metadata file exists
if (!file.exists(opt$sample_metadata_file)) {
  stop("sample metadata file not found.")
}

# read in mitochondrial gene list
mito_genes <- unique(scan(opt$mito_file, what = "character"))

# read in gtf file (genes only for speed)
gtf <- rtracklayer::import(opt$gtf_file, feature.type = "gene")

# parse sample id list
sample_ids <- unlist(stringr::str_split(opt$sample_id, ",|;")) |> sort()

# set include unspliced for non feature data
include_unspliced <- !opt$spliced_only

# get unfiltered sce
unfiltered_sce <- read_alevin(
  quant_dir = opt$alevin_dir,
  include_unspliced = include_unspliced,
  fry_mode = TRUE,
  tech_version = opt$technology,
  assay_ontology_term_id = opt$assay_ontology_term_id,
  seq_unit = opt$seq_unit,
  library_id = opt$library_id,
  sample_id = sample_ids,
  project_id = opt$project_id
)

# read and merge feature counts if present
if (opt$feature_dir != "") {
  feature_sce <- read_alevin(
    quant_dir = opt$feature_dir,
    include_unspliced = FALSE,
    fry_mode = TRUE,
    feature_data = TRUE,
    tech_version = opt$technology,
    assay_ontology_term_id = opt$assay_ontology_term_id,
    seq_unit = opt$seq_unit,
    library_id = opt$library_id,
    sample_id = sample_ids,
    project_id = opt$project_id
  )

  unfiltered_sce <- merge_altexp(unfiltered_sce, feature_sce, opt$feature_name)
  # add alt experiment features stats
  altExp(unfiltered_sce, opt$feature_name) <- scuttle::addPerFeatureQCMetrics(altExp(unfiltered_sce, opt$feature_name))

  # if CITE, add `adt_id` column to rowData with rownames
  if (opt$feature_name == "adt") {
    rowData(altExp(unfiltered_sce, "adt"))$adt_id <- rownames(rowData(altExp(unfiltered_sce, "adt")))
  }
}


# read in sample metadata
sample_metadata_df <- readr::read_tsv(opt$sample_metadata_file) |>
  # rename sample id column
  dplyr::rename("sample_id" = "scpca_sample_id") |>
  # add library ID as column in sample metadata
  # we need this so we are able to merge sample metadata with colData later
  dplyr::mutate(library_id = opt$library_id)

if ("upload_date" %in% colnames(sample_metadata_df)) {
  sample_metadata_df <- sample_metadata_df |>
    # remove upload date as we don't provide this on the portal
    dplyr::select(-upload_date)
}

# add per cell and per gene statistics to colData and rowData
unfiltered_sce <- unfiltered_sce |>
  add_cell_mito_qc(mito = mito_genes) |>
  # add gene symbols to rowData
  add_gene_symbols(gene_info = gtf) |>
  scuttle::addPerFeatureQCMetrics() |>
  # add dataframe with sample metadata to sce metadata
  # `add_sample_metadata` will filter sample_metadata_df to the relevant sample ids
  add_sample_metadata(metadata_df = sample_metadata_df)

# if columns with sample type info aren't provided, set to NA
if (!("is_xenograft" %in% colnames(sample_metadata_df))) {
  sample_metadata_df$is_xenograft <- NA
}
if (!("is_cell_line" %in% colnames(sample_metadata_df))) {
  sample_metadata_df$is_cell_line <- NA
}

# add explicit metadata field for the sample type
sample_type <- sample_metadata_df |>
  dplyr::filter(sample_id %in% sample_ids) |>
  dplyr::mutate(
    sample_type = dplyr::case_when(
      is_xenograft ~ "patient-derived xenograft",
      is_cell_line ~ "cell line",
      # if neither column was provided, note that
      is.na(is_xenograft) & is.na(is_cell_line) ~ "Not provided",
      .default = "patient tissue"
    )
  ) |>
  dplyr::select(sample_id, sample_type) |>
  # convert into named vector
  tibble::deframe()

# unname if length is 1, and add to sce metadata
if (length(sample_type) == 1) {
  sample_type <- unname(sample_type)
}
metadata(unfiltered_sce)$sample_type <- sample_type

# write to rds
readr::write_rds(unfiltered_sce, opt$unfiltered_file, compress = "bz2")
