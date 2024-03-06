#!/usr/bin/env Rscript

# import libraries
library(optparse)
library(SingleCellExperiment)

# Set up optparse options
option_list <- list(
  make_option(
    opt_str = c("--input_library_ids"),
    type = "character",
    help = "Comma separated list of library IDs corresponding to the libraries being merged."
  ),
  make_option(
    opt_str = c("--input_sce_files"),
    type = "character",
    help = "Comma separated list of input sce file paths corresponding to the sces being merged."
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
  stop("Only 1 input file provided, no merging will be performed for this group")
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

# Read in SCEs -----------------------------------------------------------------

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

# Add cell type annotation columns where needed  -------------------------------

# check for present cell type annotations
all_celltypes <- sce_list |>
  purrr::map(\(sce) {
    metadata(sce)$celltype_methods
  }) |>
  purrr::reduce(union)


if ("submitter" %in% all_celltypes) {
  # Add `"Submitter-excluded"` value to any libraries without submitter
  sce_list <- sce_list |>
    purrr::map(\(sce){
      if (!"submitter" %in% metadata(sce)$celltype_methods) {
        colData(sce)$submitter_celltype_annotation <- "Submitter-excluded"
      }
      return(sce)
    })
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

  # Add `"Cell type annotation not performed"` string to libraries without SingleR
  sce_list <- sce_list |>
    purrr::map(\(sce){
      if (!"singler" %in% metadata(sce)$celltype_methods) {
        colData(sce)$singler_celltype_annotation <- "Cell type annotation not performed"
        if (use_ontology) {
          colData(sce)$singler_celltype_ontology <- "Cell type annotation not performed"
        }
      }
      return(sce)
    })
}
if ("cellassign" %in% all_celltypes) {
  # Add `"Cell type annotation not performed"` string to libraries without CellAssign,
  #  and make the max prediction `NA_real_` for extra safety
  sce_list <- sce_list |>
    purrr::map(\(sce){
      if (!"cellassign" %in% metadata(sce)$celltype_methods) {
        colData(sce)$cellassign_celltype_annotation <- "Cell type annotation not performed"
        colData(sce)$cellassign_max_prediction <- NA_real_
      }
      return(sce)
    })
}


# Update some SCE information  -------------------------------------------------
# - Add a new colData column with any additional modalities
# - Remove cluster parameters and miQC model from metadata
sce_list <- sce_list |>
  purrr::map(\(sce){
    additional_modalities <- altExpNames(sce)
    if (length(additional_modalities) == 0) {
      additional_modalities <- NA
    }
    sce$additional_modalities <- paste0(additional_modalities, collapse = ";")

    metadata(sce)$cluster_algorithm <- NULL
    metadata(sce)$cluster_weighting <- NULL
    metadata(sce)$cluster_nn <- NULL
    metadata(sce)$miQC_model <- NULL

    return(sce)
  })

# Determine SCE columns to retain  ---------------------------------------------

# Define all possible colData columns which we might want to retain
possible_columns <- c(
  "barcodes",
  # RNA statistics and post-processing
  "sum",
  "detected",
  "total",
  "subsets_mito_sum",
  "subsets_mito_detected",
  "subsets_mito_percent",
  "miQC_pass",
  "prob_compromised",
  "scpca_filter",
  "additional_modalities",
  # cell type names
  "submitter_celltype_annotation",
  "singler_celltype_annotation",
  "singler_celltype_ontology",
  "cellassign_celltype_annotation",
  "cellassign_max_prediction",
  # ADT statistics
  "adt_scpca_filter",
  "altexps_adt_sum",
  "altexps_adt_detected",
  "altexps_adt_percent"
)

# Define colData columns to retain based on intersection with present columns
present_columns <- sce_list |>
  purrr::map(\(sce) names(colData(sce))) |>
  purrr::reduce(union)
retain_coldata_columns <- intersect(possible_columns, present_columns)

# Define altExp columns to retain/preserve, currently only for "adt"
adt_possible_columns <- c(
  "zero.ambient",
  "high.ambient",
  "ambient.scale",
  "sum.controls",
  "high.controls",
  "discard"
)
adt_present_columns <- sce_list |>
  # only consider "adt" altExps
  purrr::keep(\(sce) "adt" %in% altExpNames(sce)) |>
  purrr::map(\(sce) names(colData(altExp(sce, "adt")))) |>
  # use .init to handle an empty input; will return the .init value
  # if input is empty, then there are no "adt" altExps anyways
  purrr::reduce(union, .init = NULL)

# ensure that there are indeed no "adt" altExps if adt_present_columns is empty
adt_altexps <- sce_list |>
  purrr::keep(\(sce) "adt" %in% altExpNames(sce))
if (is.null(adt_present_columns) && length(adt_altexps) > 0) {
  stop("Error in determining which adt altExp columns should be retained.")
}

retain_altexp_coldata_list <- list("adt" = intersect(adt_possible_columns, adt_present_columns))
preserve_altexp_rowdata_list <- list("adt" = c("adt_id", "target_type"))

# Merge SCEs -------------------------------------------------------------------


# create combined SCE object
merged_sce <- scpcaTools::merge_sce_list(
  sce_list,
  batch_column = "library_id",
  cell_id_column = "cell_id",
  retain_coldata_cols = retain_coldata_columns,
  include_altexp = opt$include_altexp,
  preserve_rowdata_cols = c("gene_symbol", "gene_ids"),
  retain_altexp_coldata_cols = retain_altexp_coldata_list,
  preserve_altexp_rowdata_cols = preserve_altexp_rowdata_list
)

# add sample metadata to colData
merged_sce <- scpcaTools::metadata_to_coldata(
  merged_sce,
  join_columns = "library_id"
)

# remove sample metadata
metadata(merged_sce) <- metadata(merged_sce)[names(metadata(merged_sce)) != "sample_metadata"]

# grab technology and EFO from metadata$library_metadata
library_df <- names(sce_list) |>
  purrr::map(\(library_id){
    lib_meta <- metadata(merged_sce) |>
      purrr::pluck("library_metadata", library_id)
    data.frame(
      library_id = library_id,
      tech_version = lib_meta$tech_version,
      assay_ontology_term_id = lib_meta$assay_ontology_term_id,
      # rename to suspension_type for consistency with merged AnnData objects
      suspension_type = lib_meta$seq_unit
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

metadata(merged_sce)$merged_highly_variable_genes <- hvg_list

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
