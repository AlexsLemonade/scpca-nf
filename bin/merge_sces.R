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
  ),
  make_option(
    opt_str = c("--n_hvg"),
    type = "integer",
    default = 2000,
    help = "number of high variance genes to use for dimension reduction;
            the default is n_hvg = 2000"
  ),
  make_option(
    opt_str = c("--include_altexp"),
    action = "store_true",
    default = FALSE,
    help = "Keep any altExp present in the merged object."
  ),
  make_option(
    opt_str = c("--multiplexed"),
    action = "store_true",
    default = FALSE,
    help = "Indicates if the provided SCE's contain multiplexed data.
      If so, the sample metadata will not be added to the colData."
  ),
  make_option(
    opt_str = c("-t", "--threads"),
    type = "integer",
    default = 1,
    help = "Number of multiprocessing threads to use"
  )
)

# Setup ------------------------------------------------------------------------

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

# check that file extension for output file is correct
if (!(stringr::str_ends(opt$output_sce_file, ".rds"))) {
  stop("output file name must end in .rds")
}

# check that input files are provided and more than one are present for merging and integration
if (is.null(opt$input_sce_files)) {
  stop("List of input files containing individual SCE objects to merge is missing.")
} else {
  # list of paths to sce files
  input_sce_files <- unlist(stringr::str_split(opt$input_sce_files, ","))
}

if (length(input_sce_files) == 1) {
  stop("Only 1 input file provided, no merging or integration will be performed for this group")
}

# use library ids to name list of input files
if (is.null(opt$input_library_ids)) {
  # extract library ids from filename
  names(input_sce_files) <- stringr::word(basename(input_sce_files), 1, sep = "_")
} else {
  # pull library ids from list
  names(input_sce_files) <- unlist(stringr::str_split(opt$input_library_ids, ","))
}

# sort inputs by library id
input_sce_files <- input_sce_files[order(names(input_sce_files))]

# check that input files exist
missing_sce_files <- input_sce_files[which(!file.exists(input_sce_files))]
if (length(missing_sce_files) > 0) {
  stop(
    glue::glue(
      "\nCannot find input file: {missing_sce_files}."
    )
  )
}

# set up multiprocessing params
if (opt$threads > 1) {
  bp_param <- BiocParallel::MulticoreParam(opt$threads)
} else {
  bp_param <- BiocParallel::SerialParam()
}

# Merge SCEs -------------------------------------------------------------------

# get list of sces
sce_list <- purrr::map(input_sce_files, readr::read_rds)

# check that all input RDS files contain SCE objects
sce_checks <- purrr::map_lgl(
  sce_list,
  \(x) is(x, "SingleCellExperiment")
)
if (!all(sce_checks)) {
  stop(
    "All input files must contain a `SingleCellExperiment` object."
  )
}

# check if any SCEs have cell types; if so, we want to retain those colData columns
all_celltypes <- sce_list |>
  purrr::map(\(sce) {
    metadata(sce)$celltype_methods
  }) |>
  purrr::reduce(union)

# Default vector of colData columns to retain from
# https://github.com/AlexsLemonade/scpcaTools/blob/2ebdc4f4dfc4233fad97805f9a9a5e3bc6919f1e/R/merge_sce_list.R
retain_coldata_columns <- c(
  "sum",
  "detected",
  "total",
  "subsets_mito_sum",
  "subsets_mito_detected",
  "subsets_mito_percent",
  "miQC_pass",
  "prob_compromised",
  "barcodes"
)
# add relevant cell type columns to `retain_coldata_columns`
if ("submitter" %in% all_celltypes) {
  retain_coldata_columns <- c(
    retain_coldata_columns,
    "submitter_celltype_annotation"
  )
}
if ("singler" %in% all_celltypes) {
  # Check if the label used for annotation was ontology in at least 1 SCE
  use_ontology <- sce_list |>
    purrr::map(\(sce){
      metadata(sce)$singler_reference_label == "label.ont"
    }) |>
    # avoid warning with unlist; can't use map_lgl since ^ would always need to
    #  return length 1
    unlist() |>
    any()

  retain_coldata_columns <- c(
    retain_coldata_columns,
    "singler_celltype_annotation",
    # only use ontology if TRUE
    ifelse(use_ontology, "singler_celltype_ontology", NULL)
  )
}
if ("cellassign" %in% all_celltypes) {
  retain_coldata_columns <- c(
    retain_coldata_columns,
    "cellassign_celltype_annotation",
    "cellassign_max_prediction"
  )
}


# Add a new column with any additional modalities
sce_list <- sce_list |>
  purrr::map(\(sce){
    # value will be adt, cellhash, or NA
    additional_modalities <- altExpNames(sce)
    if (length(additional_modalities) == 0) {
      additional_modalities <- NA
    }
    sce$additional_modalities <- additional_modalities
    return(sce)
  })

# create combined SCE object
merged_sce <- scpcaTools::merge_sce_list(
  sce_list,
  batch_column = "library_id",
  retain_coldata_cols = retain_coldata_columns,
  preserve_rowdata_cols = "gene_symbol",
  cell_id_column = "cell_id",
  include_altexp = opt$include_altexp
)

# add sample metadata to colData as long as there are no multiplexed data
if (!opt$multiplexed) {
  merged_sce <- scpcaTools::metadata_to_coldata(
    merged_sce,
    join_columns = "library_id"
  )

  # remove sample metadata
  metadata(merged_sce) <- metadata(merged_sce)[names(metadata(merged_sce)) != "sample_metadata"]
}

# grab technology and EFO from metadata$library_metadata
library_df <- names(input_sce_files) |>
  purrr::map(\(library_id){
    lib_meta <- metadata(merged_sce) |>
      purrr::pluck("library_metadata", library_id)
    data.frame(
      library_id = library_id,
      tech_version = lib_meta$tech_version,
      assay_ontology_term_id = lib_meta$assay_ontology_term_id,
      seq_unit = lib_meta$seq_unit
    )
  }) |>
  dplyr::bind_rows()

# join tech and EFO with colData
colData(merged_sce) <- colData(merged_sce) |>
  as.data.frame() |>
  dplyr::left_join(library_df, by = c("library_id"), relationship = "many-to-one") |>
  DataFrame(row.names = rownames(colData(merged_sce)))

# HVG selection ----------------------------------------------------------------

# extract the column with the block variable
batch_column <- merged_sce$library_id

# model gene variance
gene_var_block <- scran::modelGeneVar(
  merged_sce,
  block = batch_column,
  BPPARAM = bp_param
)
# identify subset of variable genes
hvg_list <- scran::getTopHVGs(
  gene_var_block,
  n = opt$n_hvg
)

metadata(merged_sce)$merged_hvgs <- hvg_list

# Dim Reduction PCA and UMAP----------------------------------------------------

# multi batch PCA on merged object
multi_pca <- batchelor::multiBatchPCA(
  merged_sce,
  subset.row = hvg_list,
  batch = batch_column,
  preserve.single = TRUE,
  BPPARAM = bp_param
)

# add PCA results to merged SCE object
reducedDim(merged_sce, "PCA") <- multi_pca[[1]]

# add UMAP
merged_sce <- scater::runUMAP(
  merged_sce,
  dimred = "PCA",
  BPPARAM = bp_param
)

# write out merged sce file
readr::write_rds(
  merged_sce,
  opt$output_sce_file,
  compress = "gz",
  compression = 2
)
