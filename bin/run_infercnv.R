#!/usr/bin/env Rscript

# This script runs inferCNV on a SingleCellExperiment object and saves results to the SCE object:
# metadata field `infercnv_table`: (data frame) The full table output from the inferCNV HMM
# metadata field `infercnv_options`: (list) The options used when calling inferCNV, from the @options slot of the output object
# colData column `total_cnv`: The sum of CNV per cell, calculated from the HMM output

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
    help = "Path to the output SCE file inferCNV output"
  ),
  make_option(
    opt_str = c("--output_dir"),
    type = "character",
    default = "",
    help = "Folder to save final infercnv results"
  ),
  make_option(
    opt_str = c("--gene_order_file"),
    type = "character",
    default = "",
    help = "Path to gene order file as tab delimited .txt file with no column headers.
      Columns are: Ensembl gene id, chr, start, stop."
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
  "gene_order_file does not exist" = file.exists(opts$gene_order_file),
  "output_dir does not exist" = dir.exists(opts$output_dir),
  "output_rds was not provided" = !is.null(opts$output_rds)
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
png_file <- file.path(opts$output_dir, "infercnv.png")
infercnv_final_file <- file.path(opts$output_dir, "run.final.infercnv_obj")
metadata_file <- file.path(opts$output_dir, "map_metadata_from_infercnv.txt")

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
        out_dir = opts$output_dir, # save all intermediate files here
        denoise = TRUE, # this option is generally used in inferCNV HMM examples
        HMM = TRUE, # use the i6 HMM
        HMM_type = "i6",
        save_rds = FALSE, # don't save intermediate RDS files
        num_threads = opts$threads
      )
  },
  error = function(e) {
    # If inferCNV failed, create an empty heatmap file
    system(glue::glue("touch {png_file}"))

    # return NULL
    NULL
  }
)
# confirm png file exists
stopifnot("PNG file not created" = file.exists(png_file))

# save results SCE if inferCNV ran successfully ------------------------------
if (!is.null(infercnv_result)) {
  # confirm final infercnv object exists; the metadata table is created from this file
  stopifnot(
    "inferCNV did not write output to disc" = file.exists(infercnv_final_file)
  )

  # create wide table with barcodes and inferred CNV events
  # this will automatically create `metadata_file`
  infercnv::add_to_seurat(
    seurat_obj = NULL,
    infercnv_output_path = opts$output_dir
  )

  # save table to SCE metadata
  # note we have to read in with base R, since there are row names
  infercnv_table <- read.table(metadata_file, header = TRUE, sep = "\t") |>
    tibble::rownames_to_column(var = "barcodes")
  metadata(sce)$infercnv_table <- infercnv_table

  # save inferCNV runtime options used to SCE metadata
  metadata(sce)$infercnv_options <- infercnv_result@options

  # add a total_cnv column to colData
  total_cnv_df <- infercnv_table |>
    tidyr::pivot_longer(
      starts_with("has_cnv_"),
      names_to = "chr",
      values_to = "cnv"
    ) |>
    dplyr::group_by(barcodes) |>
    dplyr::summarize(total_cnv = sum(cnv))
  colData(sce) <- colData(sce) |>
    as.data.frame() |>
    dplyr::left_join(total_cnv_df, by = "barcodes") |>
    DataFrame(row.names = colnames(sce))
}

# export SCE -------------------
readr::write_rds(
  sce,
  opts$output_rds,
  compress = "bz2"
)


# clean up: remove all non-heatmap files from the output directory ------
remove_files <- list.files(opts$output_dir, full.names = TRUE)
remove_files <- remove_files[remove_files != png_file]
fs::file_delete(remove_files)
