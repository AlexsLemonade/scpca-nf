#!/usr/bin/env Rscript

# This script runs inferCNV on a SingleCellExperiment object and export result files:
# - inferCNV heatmap png
# - inferCNV RDS results file
# - wide metadata table with CNVs

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(optparse)
})

# set option per inferCNV:
# """Please use "options(scipen = 100)" before running infercnv if you are using
# the analysis_mode="subclusters" option or you may encounter an error while the
# hclust is being generated."""
# analysis_mode="subclusters" is the default, which we do use
options(scipen = 100)

option_list <- list(
  make_option(
    opt_str = "--input_sce_file",
    type = "character",
    default = "",
    help = "Path to the SCE file to run inferCNV on"
  ),
  make_option(
    opt_str = "--output_rds",
    type = "character",
    help = "Path to the output RDS file to hold the final inferCNV result objects"
  ),
  make_option(
    opt_str = "--output_table",
    type = "character",
    help = "Path to the output TSV file to hold the inferCNV metadata table of CNVs"
  ),
  make_option(
    opt_str = "--output_heatmap",
    type = "character",
    help = "Path to the PNG file to hold the inferCNV heatmap"
  ),
  make_option(
    opt_str = c("--gene_order_file"),
    type = "character",
    default = "",
    help = "Path to gene order file as tab delimited .txt file with no column headers.
      Columns are: Ensembl gene id, chr, start, stop."
  ),
  make_option(
    opt_str = c("--temp_dir"),
    type = "character",
    default = "infercnv_tmp",
    help = "Temporary directory to save inferCNV output to"
  ),
  make_option(
    opt_str = c("-t", "--threads"),
    type = "integer",
    default = 4,
    help = "Number of multiprocessing threads to use."
  ),
  make_option(
    opt_str = c("--random_seed"),
    type = "integer",
    default = 2025,
    help = "Random seed to set for reproducibility. Note that inferCNV is only reproducible on a given operating system."
  )
)

# parse and check input options ----------------
opts <- parse_args(OptionParser(option_list = option_list))

stopifnot(
  "input_sce_file does not exist" = file.exists(opts$input_sce_file),
  "output_rds was not provided" = !is.null(opts$output_rds),
  "output_table was not provided" = !is.null(opts$output_table),
  "output_heatmap was not provided" = !is.null(opts$output_heatmap),
  "gene_order_file does not exist" = file.exists(opts$gene_order_file)
)

# set up files and paths -----------------
set.seed(opts$seed)

# read SCE file and check for column
sce <- readRDS(opts$input_sce_file)
stopifnot(
  "SCE colData is missing `is_infercnv_reference` column" =
    "is_infercnv_reference" %in% colnames(colData(sce))
)

# define relevant infercnv output files for later use/checks
# infercnv will automatically create these files at these hardcoded paths
fs::dir_create(opts$temp_dir) # ensure output directory exists, to be safe
scratch_infercnv_rds <- file.path(opts$temp_dir, "run.final.infercnv_obj")
scratch_metadata_file <- file.path(opts$temp_dir, "map_metadata_from_infercnv.txt")
scratch_png_file <- file.path(opts$temp_dir, "infercnv.png")

# run infercnv ------------------------

# create annotations_df table which requires:
# - rownames as cell barcodes
# - a single column with annotation labels; we use reference and query
# see `data("infercnv_annots_example")` for context
annotations_df <- data.frame(
  annotation = ifelse(sce$is_infercnv_reference, "reference", "query"),
  row.names = sce$barcodes
)

# run infercnv in a tryCatch in case of errors
# this will automatically create the heatmap file in the output directory
infercnv_result <- tryCatch(
  {
    # create the inferCNV object and pipe into run
    # so errors from either step are caught
    infercnv::CreateInfercnvObject(
      raw_counts_matrix = counts(sce),
      annotations_file = annotations_df,
      gene_order_file = opts$gene_order_file,
      ref_group_name = "reference",
      # use our chr names with arms; chrM was already removed
      chr_exclude = c("chrXp", "chrXq", "chrYp", "chrYq")
    ) |>
      infercnv::run(
        cutoff = 0.1, # use 1 for smart-seq, 0.1 for 10x-genomics
        out_dir = opts$temp_dir, # save all intermediate files here
        denoise = TRUE, # this option is generally used in inferCNV HMM examples
        HMM = TRUE, # use the i6 HMM
        HMM_type = "i6",
        save_rds = FALSE, # don't save intermediate RDS files
        num_threads = opts$threads
      )
  },
  error = function(e) {
    message("inferCNV failed; creating empty result files")

    # If inferCNV failed, create empty result files
    file.create(
      opts$output_rds,
      opts$output_table,
      opts$output_heatmap
    )

    # return NULL
    NULL
  }
)

# save relevant results to RDS if inferCNV ran successfully -------------------
if (!is.null(infercnv_result)) {
  # confirm final infercnv object exists
  stopifnot(
    "inferCNV did not write expected output" = file.exists(scratch_infercnv_rds)
  )

  # create wide table with barcodes and inferred CNV events
  # this will automatically create `scratch_metadata_file`
  infercnv::add_to_seurat(
    seurat_obj = NULL,
    infercnv_output_path = opts$temp_dir
  )

  # rename result files to final names
  fs::file_move(scratch_infercnv_rds, opts$output_rds)
  fs::file_move(scratch_metadata_file, opts$output_table)
  fs::file_move(scratch_png_file, opts$output_heatmap)
}


# confirm all final files exist
stopifnot(
  "inferCNV results file not created" = file.exists(opts$output_rds),
  "inferCNV metadata table file not created" = file.exists(opts$output_table),
  "PNG file not created" = file.exists(opts$output_heatmap)
)
