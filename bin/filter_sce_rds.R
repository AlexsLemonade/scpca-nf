#!/usr/bin/env Rscript

# This script takes a SingleCellExperiment stored in a .rds file and
# filters it using emptyDrops, adding miQC metrics for probability compromised

# import libraries
suppressPackageStartupMessages({
  library(optparse)
  library(SingleCellExperiment)
})
# set up arguments
option_list <- list(
  make_option(
    opt_str = c("-u", "--unfiltered_file"),
    type = "character",
    help = "path to rds file with unfiltered sce object"
  ),
  make_option(
    opt_str = c("-f", "--filtered_file"),
    type = "character",
    help = "path to output filtered rds file. Must end in .rds"
  ),
  make_option(
    opt_str = c("--prob_compromised_cutoff"),
    type = "double",
    default = 0.75,
    help = "probability compromised cutoff used for filtering cells with miQC"
  ),
  make_option(
    opt_str = c("--enforce_left_cutoff"),
    action = "store_true",
    default = FALSE,
    help = "flag to set the miQC's enforce_left_cutoff option to TRUE"
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
   opt_str = c("-r", "--random_seed"),
   type = "integer",
   help = "A random seed for reproducibility."
 )
)

opt <- parse_args(OptionParser(option_list = option_list))

# set seed
set.seed(opt$random_seed)

# check that unfiltered file file exists
if(!file.exists(opt$unfiltered_file)){
  stop("Missing unfiltered.rds file")
}

# check that output file name ends in .rds
if(!(stringr::str_ends(opt$filtered_file, ".rds"))){
  stop("filtered file name must end in .rds")
}

# check that prob compromised cutoff is between 0-1
if(!dplyr::between(opt$prob_compromised_cutoff, 0, 1)){
  stop("--prob_compromised_cutoff must be a number between 0 to 1")
}

# read in unfiltered rds file
unfiltered_sce <- readr::read_rds(opt$unfiltered_file)

# check if ADT data is present, and prep if so
if (is.null(opt$adt_barcode_file)) {
  ambient_profile <- NULL
} else {
  # assign and check name for this alternative experiment
  adt_exp <- opt$adt_name
  if (!adt_exp %in% altExpNames(unfiltered_sce)) {
    stop("Given named ADT alternative experiment not present in unfiltered SCE.")
  }

  # Calculate ambient profile from empty drops for later use, and save to altExp metadata
  ambient_profile <- DropletUtils::ambientProfileEmpty( counts(altExp(unfiltered_sce, adt_exp)) )
  metadata(altExp(unfiltered_sce, adt_exp))$ambient_profile <- ambient_profile
}

# filter sce
filtered_sce <- scpcaTools::filter_counts(unfiltered_sce)
# remove unfiltered for memory saving
rm(unfiltered_sce)

# need to remove old gene-level rowData statistics first
drop_cols <- colnames(rowData(filtered_sce)) %in% c('mean', 'detected')
rowData(filtered_sce) <- rowData(filtered_sce)[!drop_cols]

# recalculate rowData and add to filtered sce
filtered_sce <- filtered_sce |>
  scuttle::addPerFeatureQCMetrics()

# add prob_compromised to colData and miQC model to metadata
# since this can fail, we will check for success
miQC_worked <- FALSE
try({
  filtered_sce <- scpcaTools::add_miQC(filtered_sce,
                                       posterior_cutoff = opt$prob_compromised_cutoff,
                                       enforce_left_cutoff = opt$enforce_left_cutoff)
  metadata(filtered_sce)$prob_compromised_cutoff <- opt$prob_compromised_cutoff
  miQC_worked <- TRUE
 })

# set prob_compromised to NA if miQC failed
if (!miQC_worked){
  warning("miQC failed. Setting `prob_compromised` to NA.")
  filtered_sce$prob_compromised <- NA_real_
  filtered_sce$miQC_pass <- NA
  metadata(filtered_sce)$prob_compromised_cutoff <- NA
}

# grab names of altExp, if any
alt_names <- altExpNames(filtered_sce)

for (alt in alt_names) {
  # remove old row data from unfiltered
  drop_cols <- colnames(rowData(altExp(filtered_sce, alt))) %in% c('mean', 'detected')
  rowData(altExp(filtered_sce, alt)) <- rowData(altExp(filtered_sce, alt))[!drop_cols]

  # add alt experiment features stats for filtered data
  altExp(filtered_sce, alt) <- scuttle::addPerFeatureQCMetrics(altExp(filtered_sce, alt))
}


# calculate filtering QC from ADTs, if present
if (!is.null(ambient_profile)) {

  # Create data frame of ADTs and their target types
  # If `target_type` column is not present, assume all ADTs are targets
  adt_barcode_df <- readr::read_tsv(
    opt$adt_barcode_file,
    # if only 2 columns exist, only the first two col_names will be used
    col_names = c("name", "barcode", "target_type")
  )

  # Assign default of `target` if no third column provided
  if (!"target_type" %in% names(adt_barcode_df)) {
    adt_barcode_df$target_type <- "target"
  }
  # Check for valid target types and warn if invalid found
  allowed_adt_values <- c("target", "neg_control", "pos_control")
  if (!all(adt_barcode_df$target_type) %in% allowed_adt_values) {
    warn("Warning: Invalid target type values provided in the ADT barcode file.
         These ADTs will be assumed to be targets.")

    # Set all invalid values to "target" default
    adt_barcode_df <- adt_barcode_df |>
      dplyr::mutate(target_type == ifelse(
        target_type %in% allowed_adt_values,
        target_type,
        "target"
      ))
  }

  # Check that barcode file ADTs match SCE ADTs
  if ( !all.equal( sort(adt_barcode_df$name), sort(rownames(altExp(filtered_sce, adt_exp))) )) {
    stop("Mismatch between provided ADT barcode file and ADTs in SCE.")
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

  # Save QC stats to the altexp
  colData(altExp(filtered_sce, adt_exp)) <- cbind(colData(altExp(filtered_sce, adt_exp)),
                                                  adt_qc_df)

  # Add `target_type` to rowData
  rowData(altExp(filtered_sce, adt_exp))$target_type <- adt_barcode_df$target_type
}


# write filtered sce to output
readr::write_rds(filtered_sce, opt$filtered_file, compress = "gz")

