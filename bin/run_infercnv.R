#!/usr/bin/env Rscript

# This script runs inferCNV on a SingleCellExperiment object and export result files:
# - inferCNV heatmap png
# - an RDS file with a list of information to save to the SCE file later:
#  - wide metadata table with CNVs
#  - list of inferCNV options

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(optparse)
})

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
    help = "Path to the output RDS file to hold inferCNV results of interest"
  ),
  make_option(
    opt_str = "--output_heatmap",
    type = "character",
    help = "Path to the output heatmap PNG file"
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
  "output_heatmap was not provided" = !is.null(opts$output_heatmap),
  "temp_dir was not provided" = !is.null(opts$temp_dir),
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
scratch_png_file <- file.path(opts$temp_dir, "infercnv.png")
scratch_infercnv_final_file <- file.path(opts$temp_dir, "run.final.infercnv_obj")
scratch_metadata_file <- file.path(opts$temp_dir, "map_metadata_from_infercnv.txt")

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
    message("inferCNV failed; creating an empty heatmap")

    # If inferCNV failed, create an empty heatmap file at the _final destination_
    system(glue::glue("touch {opts$output_heatmap}"))

    # return NULL
    NULL
  }
)

# save relevant results to RDS if inferCNV ran successfully -------------------
if (!is.null(infercnv_result)) {
  # confirm final infercnv object exists; the metadata table is created from this file
  stopifnot(
    "inferCNV did not write output to disc" = file.exists(scratch_infercnv_final_file)
  )

  # create wide table with barcodes and inferred CNV events
  # this will automatically create `scratch_metadata_file`
  infercnv::add_to_seurat(
    seurat_obj = NULL,
    infercnv_output_path = opts$temp_dir
  )

  # create list of results to export

  # note we have to read in with base R, since there are row names
  infercnv_table <- read.table(scratch_metadata_file, header = TRUE, sep = "\t") |>
    tibble::rownames_to_column(var = "barcodes")

  output_results <- list(
    infercnv_table = infercnv_table,
    infercnv_options = infercnv_result@options
  )

  # export results
  readr::write_rds(output_results, opts$output_rds, compress = "bz2")

  # copy heatmap to final destination
  fs::file_copy(scratch_png_file, opts$output_heatmap, overwrite = TRUE)
}


# confirm final png file exists
stopifnot("PNG file not created" = file.exists(opts$output_heatmap))
