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
    opt_str = c("-i", "--filtered_sce_file"),
    type = "character",
    help = "path to rds file with input sce object to be processed"
  ),
  make_option(
    opt_str = c("-o", "--output_sce_file"),
    type = "character",
    help = "path to output rds file to store processed sce object. Must end in .rds"
  ),
  make_option(
    opt_str = c("--adt_name"),
    type = "character",
    default = "adt",
    help = "Name for the alternative experiment, if present, that contains ADT features"
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

# check that filtered SCE file exists
if(!file.exists(opt$filtered_sce_file)){
  stop("Missing filtered.rds file")
}

# check that output file name ends in .rds
if(!(stringr::str_ends(opt$output_sce_file, ".rds"))){
  stop("output file name must end in .rds")
}

# read in filtered rds file
sce <- readr::read_rds(opt$filtered_sce_file)

# create scpca_filter column -----------
if(all(is.na(sce$miQC_pass))){
  sce$scpca_filter <- ifelse(
    sce$detected >= opt$gene_cutoff,
    "Keep",
    "Remove"
  )
  metadata(sce)$scpca_filter_method <- "Minimum_gene_cutoff"
} else {
  # create filter using miQC and min gene cutoff
  sce$scpca_filter <- ifelse(
    sce$miQC_pass & sce$detected >= opt$gene_cutoff,
    "Keep",
    "Remove"
  )
  metadata(sce)$scpca_filter_method <- "miQC"
}

# add min gene cutoff to metadata
metadata(sce)$min_gene_cutoff <- opt$gene_cutoff

# create adt_scpca_filter column, if CITEseq ----------
alt_exp <- opt$adt_name
if (alt_exp %in% altExpNames(sce)) {
  
  # check for NAs in filtering columns
  use_discard <- sum(is.na(altExp(sce, alt_exp)$discard)) == 0
  use_zero.ambient <- sum(is.na(altExp(sce, alt_exp)$zero.ambient)) == 0
  
  # Fail right away if both are FALSE
  if (!use_discard & !use_zero.ambient) {
    warning("Failed to filter ADT data.")
    sce$adt_scpca_filter <- NULL
    metadata(sce)$adt_scpca_filter_method <- "No filter"
  } else {
    # Set up filtering based on either discard or zero.ambient
    if (use_discard) {
      sce$adt_scpca_filter <-  ifelse(
        altExp(sce, alt_exp)$discard,
        "Remove",
        "Keep"
      )
    } else if (use_zero.ambient) {
      # use _opposite_ logic from `discard` column
      sce$adt_scpca_filter <-  ifelse(
        altExp(sce, alt_exp)$zero.ambient,
        "Keep",
        "Remove"
      )
    }
    # Now that filtering has actually been applied, assign 
    #  `adt_scpca_filter_method` metadata based on colData contents
    if ("sum.controls" %in% names(colData(altExp(sce, alt_exp)))) {
      metadata(sce)$adt_scpca_filter_method <- "cleanTagCounts with isotype controls"
    } else if ("ambient.scale" %in% names(colData(altExp(sce, alt_exp)))) {
      metadata(sce)$adt_scpca_filter_method <- "cleanTagCounts without isotype controls"
    } else {
      stop("Error in ADT filtering.")
    }
  }

}

# Perform filtering -----------------

# filter sce using criteria in scpca_filter (not adt_scpca_filter)
processed_sce <- sce[, which(sce$scpca_filter == "Keep")]

# replace existing stats with recalculated gene stats
drop_cols <- colnames(rowData(processed_sce, alt)) %in% c('mean', 'detected')
rowData(processed_sce) <- rowData(processed_sce)[!drop_cols]

processed_sce <- processed_sce |>
  scuttle::addPerFeatureQCMetrics()

# replace existing stats from altExp if any
for (alt in altExpNames(processed_sce)) {
  # remove old row data
  drop_cols <- colnames(rowData(altExp(processed_sce, alt))) %in% c('mean', 'detected')
  rowData(altExp(processed_sce, alt)) <- rowData(altExp(processed_sce, alt))[!drop_cols]

  # add alt experiment features stats for filtered data
  altExp(processed_sce, alt) <- scuttle::addPerFeatureQCMetrics(altExp(processed_sce, alt))
}

# Perform normalization -----------------

# cluster prior to RNA normalization
qclust <- NULL
try({
  # try and cluster similar cells
  # clustering will fail if < 100 cells in dataset
  qclust <- scran::quickCluster(processed_sce)

  # Compute sum factors for each cell cluster grouping
  processed_sce <- scran::computeSumFactors(processed_sce, clusters = qclust)

  # Include note in metadata re: clustering before computing sum factors
  metadata(processed_sce)$normalization <- "deconvolution"
})

if (is.null(qclust)) {
  # Include note in metadata re: failed clustering
  metadata(processed_sce)$normalization <- "log-normalization"
}

# Normalize and log transform RNA counts
processed_sce <- scuttle::logNormCounts(processed_sce)

# Try to normalize ADT counts, if present
if (alt_exp %in% altExpNames(processed_sce)) {
  
  # need `all()` since,if present, this is an array
  if( !all(is.null(metadata(altExp(processed_sce, alt_exp))$ambient_profile))){
    # Calculate median size factors from the ambient profile
    altExp(processed_sce, alt_exp) <- scuttle::computeMedianFactors(
      altExp(processed_sce, alt_exp),
      reference = metadata(altExp(processed_sce, alt_exp))$ambient_profile
    )
  } else {
    # if ambient profile is not present, set sizeFactor to 0 for later warning.
    altExp(processed_sce, alt_exp)$sizeFactor <- 0
  }
  
  adt_sce <- altExp(processed_sce, alt_exp)
  # Perform filtering specifically to allow for normalization
  if (!is.null(processed_sce$adt_scpca_filter)) {
    adt_sce <- adt_sce[,processed_sce$adt_scpca_filter == "Keep"]
  }
  
  # If any size factors are not positive or filtering is null, simply use log1p
  if ( any( adt_sce$sizeFactor <= 0 )  | is.null(processed_sce$adt_scpca_filter)) {
    metadata(processed_sce)$adt_normalization <- "log-normalization"
    logcounts(adt_sce) <- log1p(counts(adt_sce))
  } else {
    # Apply normalization using size factors
    metadata(processed_sce)$adt_normalization <- "median-based"
    adt_sce <- scuttle::logNormCounts(adt_sce)    
  }
  # Now that we have logcounts, add back to `processed_sce` but
  #   with NA values for cells not included in normalization
  
  # first, get the counts matrix and make it NA
  result_matrix <- counts(altExp(processed_sce, alt_exp))
  result_matrix[,] <- NA
  
  # now get the computed logcounts & fill them in
  result_matrix[, colnames(adt_sce)] <- logcounts(adt_sce)
  
  # Check correct number of NAs:
  observed_na_count <- sum(is.na(result_matrix))
  expected_na_count <- nrow(adt_sce) * (ncol(altExp(processed_sce, alt_exp)) - ncol(adt_sce))
  if (observed_na_count < expected_na_count) {
    stop("Incorrect number of NAs recovered during ADT normalization.")
  }

  # Add result_matrix back into correct SCE as logcounts assay
  logcounts(altExp(processed_sce, alt_exp)) <- result_matrix
}


# Perform dimension reduction --------------------

# model gene variance using `scran:modelGeneVar()`
gene_variance <- scran::modelGeneVar(processed_sce)

# select the most variable genes
var_genes <- scran::getTopHVGs(gene_variance, n = opt$n_hvg)

# save the most variable genes to the metadata
metadata(processed_sce)$highly_variable_genes <- var_genes

# dimensionality reduction
# highly variable genes are used as input to PCA
processed_sce <- scater::runPCA(processed_sce,
                               ncomponents = opt$n_pcs,
                               subset_row = var_genes)

# calculate a UMAP matrix using the PCA results
try({
  processed_sce <- scater::runUMAP(processed_sce,
                                  dimred = "PCA")
})

# Export --------------

# write out _original_ filtered SCE with additional filtering column
readr::write_rds(sce, opt$filtered_sce_file, compress = "gz")

# write out processed SCE
readr::write_rds(processed_sce, opt$output_sce_file, compress = "gz")

