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
    opt_str = c("-u", "--unfiltered_sce_file"),
    type = "character",
    default = NULL,
    help = "path to rds file with fully unfiltered (contains empty drops) input sce object. 
      This argument is only used if there is CITEseq data."
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
  make_option(
    opt_str = c("--adt_name"),
    type = "character",
    default = "CITEseq",
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


# check for ADT files if the barcode file is given
# if all present, prepare ambient profile
if (is.null(opt$adt_barcode_file)) {
  ambient_profile <- NULL
} else { 
  if (!file.exists(opt$adt_barcode_file)) {
    stop("adt_barcode_file does not exist")
  }
  if (!file.exists(opt$unfiltered_sce_file)) {
    stop("Missing unfiltered.rds file")
  }
  
  # Create data frame of ADTs and their target types
  # If `target_type` column is not present, assume all ADTs are targets
  adt_barcode_df <- readr::read_tsv(
    opt$adt_barcode_file, 
    # if only 2 columns exist, only the first two col_names will be used
    col_names = c("name", "barcode", "target_type")
  ) 
  if (!"target_type" %in% names(adt_barcode_df)) {
    adt_barcode_df$target_type <- "target"
  } 
  
  # Name for this alternative experiment
  adt_exp <- opt$adt_name
  
  # Read in unfiltered sce, check for adt_name
  sce <- readr::read_rds(opt$unfiltered_sce_file)
  if (!adt_exp %in% altExpNames(sce)) {
    stop("Given named ADT alternative experiment not present in unfiltered SCE.")
  }
  
  # Calculate ambient profile from empty drops for later use
  ambient_profile <- DropletUtils::ambientProfileEmpty( counts(altExp(sce, adt_exp)) )
  rm(sce)
}

# read in filtered rds file
sce <- readr::read_rds(opt$filtered_sce_file)

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

# Perform filtering -----------------

# filter sce using criteria in scpca_filter
filtered_sce <- sce[, which(sce$scpca_filter == "Keep")]

# replace existing stats with recalculated gene stats
drop_cols <- colnames(rowData(filtered_sce, alt)) %in% c('mean', 'detected')
rowData(filtered_sce) <- rowData(filtered_sce)[!drop_cols]

filtered_sce <- filtered_sce |>
  scuttle::addPerFeatureQCMetrics()

# replace existing stats from altExp if any
for (alt in altExpNames(filtered_sce)) {
  # remove old row data
  drop_cols <- colnames(rowData(altExp(filtered_sce, alt))) %in% c('mean', 'detected')
  rowData(altExp(filtered_sce, alt)) <- rowData(altExp(filtered_sce, alt))[!drop_cols]

  # add alt experiment features stats for filtered data
  altExp(filtered_sce, alt) <- scuttle::addPerFeatureQCMetrics(altExp(filtered_sce, alt))
}

# filter sce based on ADTs, if present
if (!is.null(ambient_profile)) {
  
  # Again, check for ADT name in this sce
  if (!adt_exp %in% altExpNames(filtered_sce)) {
    stop("Given named ADT alternative experiment not present in filtered SCE.")
  }
  
  # Find any negative controls
  neg_controls <- adt_barcode_df |>
    dplyr::filter(target_type == "neg_control") |>
    dplyr::pull(name)

  # Calculate QC stats, providing negative controls if present
  # note: function fails if controls is length 0 or null, so keep the `if`
  if (length(neg_controls) == 0) {
    adt_qc_df <- DropletUtils::cleanTagCounts(
      counts(altExp(filtered_sce, adt_exp)),
      ambient = ambient_profile
    )
  } else {
    adt_qc_df <- DropletUtils::cleanTagCounts(
      counts(altExp(filtered_sce, adt_exp)),
      ambient = ambient_profile, 
      controls = neg_controls
    )
  }

  # Save QC stats to the altexp sce
  colData(altExp(filtered_sce, adt_exp)) <- cbind(colData(altExp(filtered_sce, adt_exp)), adt_qc_df)
  
  # Filter according to `discard`
  filtered_sce <- filtered_sce[, which(!(adt_qc_df$discard))]

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

# write out _original_ filtered SCE with additional filtering column
readr::write_rds(sce, opt$filtered_sce_file, compress = "gz")

# write out processed SCE
readr::write_rds(filtered_sce, opt$output_sce_file, compress = "gz")

