#!/usr/bin/env Rscript

# This script takes an unfiltered SingleCellExperiment object and adds additional information
# sample metadata is added to metadata
# cell qc is added to colData

# import libraries
suppressPackageStartupMessages({
  library(optparse)
  library(SingleCellExperiment)
  library(scpcaTools)
})
# set up arguments
option_list <- list(
  make_option(
    opt_str = c("--sce_file"),
    type = "character",
    default = "",
    help = "path to unfiltered sce file to modify"
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
    opt_str = c("--sample_metadata_file"),
    type = "character",
    help = "path to tsv file containing sample metadata"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

# check inputs and outputs
stopifnot(
  "input SCE file not found" = file.exists(opt$input_sce_file),
  "Mitochondrial gene list file not found" = file.exists(opt$mito_file),
  "gtf file not found" = file.exists(opt$gtf_file),
  "sample metadata file not found" = file.exists(opt$sample_metadata_file)
)

# read in sce
unfiltered_sce <- readr::read_rds(opt$input_sce_file)

# read in mitochondrial gene list
mito_genes <- unique(scan(opt$mito_file, what = "character"))

# read in gtf file (genes only for speed)
gtf <- rtracklayer::import(opt$gtf_file, feature.type = "gene")

# parse sample id list
sample_ids <- unlist(stringr::str_split(opt$sample_id, ",|;")) |> sort()

# read in sample metadata
sample_metadata_df <- readr::read_tsv(opt$sample_metadata_file, col_types = "c") |>
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
      # account for sample being both cell line and PDX
      is_xenograft & is_cell_line ~ "cell line;patient-derived xenograft",
      is_xenograft ~ "patient-derived xenograft",
      is_cell_line ~ "cell line",
      # if neither column was provided, note that
      is.na(is_xenograft) & is.na(is_cell_line) ~ "Not provided",
      .default = "patient tissue"
    )
  ) |>
  dplyr::select(sample_id, sample_type) |>
  # convert into named list
  tibble::deframe() |>
  strsplit(split = ";")

# unname if length is 1, and add to sce metadata
if (length(sample_type) == 1) {
  sample_type <- unlist(sample_type) |>
    unname()
}
metadata(unfiltered_sce)$sample_type <- sample_type

# write to rds
readr::write_rds(unfiltered_sce, opt$input_sce_file, compress = "bz2")
