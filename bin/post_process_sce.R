#!/usr/bin/env Rscript

# This script takes a SingleCellExperiment stored in a .rds file with empty droplets removed
# and removes low quality cells, performs normalization, and dimensionality reduction.

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
    help = "path to rds file with input sce object to be processed"
  ),
  make_option(
    opt_str = c("-o", "--output_sce_file"),
    type = "character",
    help = "path to output rds file to store processed sce object. Must end in .rds"
  ),
  make_option(
    opt_str = c("--adt_barcode_file"),
    type = "character",
    default = NULL,
    help = "Optional path to an ADT barcode file, where the third column indicates
      the ADT type (target, negative control, or positive control)."
  ),
  # TODO: keep or remove?
  make_option(
    opt_str = c("--adt_name"),
    type = "character",
    default = "ALT", # same default as used in `bin/generate_unfiltered_sce.R`
    help = "Name for the alternative experiment, if present, that contains ADT features"
  ),
  make_option(
    opt_str = c("--adt_negative_cutoff"),
    type = "integer",
    default = 25, # TODO: pick a number..?
    help = "Maximum count of negative control antibodies to use when filtering ADT counts, if present."
  ),
  make_option(
    opt_str = c("--gene_cutoff"),
    type = "integer",
    default = 200,
    help = "Minimum number of genes per cell cutoff used for filtering cells."
  ),
  make_option(
    opt_str = c("--n_hvg"),
    type = "integer",
    default = 2000,
    help = "number of high variance genes to use for dimension reduction;
            the default is n_hvg = 2000",
  ),
  make_option(
    opt_str = c("--n_pcs"),
    type = "integer",
    default = 50,
    help = "Number of principal components to retain in the returned SingleCellExperiment object
      when using scater::runPCA; the default is n_pcs = 50",
  ),
  make_option(
    opt_str = c("-r", "--random_seed"),
    type = "integer",
    help = "A random seed for reproducibility."
  )
)

# Parse and apply input options -------------
opt <- parse_args(OptionParser(option_list = option_list))

# set seed
set.seed(opt$random_seed)

# check that input file file exists
if(!file.exists(opt$input_sce_file)){
  stop("Missing unfiltered.rds file")
}

# check that output file name ends in .rds
if(!(stringr::str_ends(opt$output_sce_file, ".rds"))){
  stop("output file name must end in .rds")
}

# read in input rds file
sce <- readr::read_rds(opt$input_sce_file)

# create scpca_filter column
if(all(is.na(sce$miQC_pass))){
  colData(sce)$scpca_filter <- ifelse(
    sce$detected >= opt$gene_cutoff, "Keep", "Remove"
  )
  metadata(sce)$scpca_filter_method <- "Minimum_gene_cutoff"
} else {
  # create filter using miQC and min gene cutoff
  colData(sce)$scpca_filter <- ifelse(
    sce$miQC_pass & sce$detected >= opt$gene_cutoff,
    "Keep",
    "Remove"
  )
  metadata(sce)$scpca_filter_method <- "miQC"
}

# add min gene cutoff to metadata
metadata(sce)$min_gene_cutoff <- opt$gene_cutoff

# read in ADT feature barcode file, if present
if (!(is.null(opt$adt_barcode_file))) {

  # Create data frame of ADTs and their target types
  adt_barcode_df <- readr::read_tsv(opt$adt_barcode_file, col_names=FALSE)

  # Check if 3rd column is present and assign column names, values assuming
  #   that all ADTs are targets if there is no 3rd column
  if (ncol(adt_barcode_df) == 2) {
    colnames(adt_barcode_df) <- c("name", "barcode")
    adt_barcode_df$target_type <- "target"
  } else if (ncol(adt_barcode_df) == 3) {
    colnames(adt_barcode_df) <- c("name", "barcode", "target_type")
  } else {
    stop("The ADT feature barcode file must have either 2 or 3 columns.")
  }

  # Name for this alternative experiment.
  altexp_name <- opt$adt_name

  # Calculate ambient profile from unfiltered sce for later use
  ambient_profile <- DropletUtils::ambientProfileEmpty( counts(altExp(sce, altexp_name)) )

  # define indicator for if ADTs need to be processed
  process_adt <- TRUE
} else {
  process_adt <- FALSE
}

# Perform filtering -----------------

# filter sce using criteria in scpca_filter
filtered_sce <- sce[, which(sce$scpca_filter == "Keep")]

filtered_sce <- filtered_sce |>
  scuttle::addPerFeatureQCMetrics()

# replace existing stats with recalculated gene stats
drop_cols = colnames(rowData(filtered_sce, alt)) %in% c('mean', 'detected')
rowData(filtered_sce) <- rowData(filtered_sce)[!drop_cols]

# replace existing stats from altExp if any
for (alt in altExpNames(filtered_sce)) {
  # remove old row data
  drop_cols <- colnames(rowData(altExp(filtered_sce, alt))) %in% c('mean', 'detected')
  rowData(altExp(filtered_sce, alt)) <- rowData(altExp(filtered_sce, alt))[!drop_cols]

  # add alt experiment features stats for filtered data
  altExp(filtered_sce, alt) <- scuttle::addPerFeatureQCMetrics(altExp(filtered_sce, alt))
}

# filter sce based on ADTs
if (process_adt) {

  adt_qc_df <- DropletUtils::cleanTagCounts(
    counts(altExp(filtered_sce, altexp_name)),
    ambient = ambient_profile
  )

  # Filter according to `discard`
  filtered_sce <- filtered_sce[, which(!(adt_qc_df$discard))]

  # TODO: Should we add `adt_qc_df` into the altExp? If we want it:
  #colData(altExp(filtered_sce, altexp_name)) <- cbind(colData(altExp(filtered_sce, altexp_name)), adt_qc_df)

  # Filter on negative controls, if present
  neg_controls <- adt_barcode_df |>
    dplyr::filter(target_type == "neg_control") |>
    dplyr::pull(name)

  if (length(neg_controls) > 0) {
    # counts for each negative control
    neg_control_counts <- counts(altExp(filtered_sce, altexp_name))[neg_controls, ]

    # remove any cell with >=threshold for any given negative control
    # first, determine which cells need to be removed
    cells_to_remove <- purrr::map(
      neg_controls,
      \(x) which(neg_control_counts[x,] >= opt$adt_negative_cutoff)
    ) |>
    unlist() |>
    unique()

    filtered_sce <- filtered_sce[, -cells_to_remove]
  }

}


# Perform normalization -----------------

# cluster prior to normalization
qclust <- NULL
try({
  # try and cluster similar cells
  # clustering will fail if < 100 cells in dataset
  qclust <- scran::quickCluster(filtered_sce)

  # Compute sum factors for each cell cluster grouping
  filtered_sce <- scran::computeSumFactors(filtered_sce, clusters = qclust)

  # Include note in metadata re: clustering before computing sum factors
  metadata(filtered_sce)$normalization <- "deconvolution"
})

if (is.null(qclust)) {
  # Include note in metadata re: failed clustering
  metadata(filtered_sce)$normalization <- "log-normalization"
}

# Normalize and log transform
filtered_sce <- scater::logNormCounts(filtered_sce)


# Perform dimension reduction --------------------

# model gene variance using `scran:modelGeneVar()`
gene_variance <- scran::modelGeneVar(filtered_sce)

# select the most variable genes
var_genes <- scran::getTopHVGs(gene_variance, n = opt$n_hvg)

# save the most variable genes to the metadata
metadata(filtered_sce)$highly_variable_genes <- var_genes

# dimensionality reduction
# highly variable genes are used as input to PCA
filtered_sce <- scater::runPCA(filtered_sce,
                               ncomponents = opt$n_pcs,
                               subset_row = var_genes)

# calculate a UMAP matrix using the PCA results
try({
  filtered_sce <- scater::runUMAP(filtered_sce,
                                  dimred = "PCA")
})

# Export --------------

# write out original SCE with additional filtering column
readr::write_rds(sce, opt$input_sce_file, compress = "gz")

# write out processed SCE
readr::write_rds(filtered_sce, opt$output_sce_file, compress = "gz")

