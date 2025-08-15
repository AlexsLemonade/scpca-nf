#!/usr/bin/env Rscript

# This script prepares and export a gene order file with chromosome arm boundaries for use with inferCNV
# This script was adapted from:
# https://github.com/AlexsLemonade/OpenScPCA-nf/blob/d314a9fc9aa8f52e938755c03a5303f0a9e48eb0/modules/infercnv-gene-order/resources/usr/bin/prepare-gene-order-files.R


library(optparse)

option_list <- list(
  make_option(
    opt_str = c("--gtf_file"),
    type = "character",
    default = "~/Desktop/Homo_sapiens.GRCh38.104.gtf.gz",
    help = "Path to input GTF file"
  ),
  make_option(
    opt_str = c("--cytoband_file"),
    type = "character",
    default = "~/Desktop/Homo_sapiens.GRCh38.104_cytoband.txt.gz",
    help = "Path to input cytoband file"
  ),
  make_option(
    opt_str = c("--gene_order_file"),
    type = "character",
    default = "homo.txt.gz",
    help = "Output file name for gene order file, which should end in `.txt.gz`."
  )
)

# Parse options ----------------------
opts <- parse_args(OptionParser(option_list = option_list))
stopifnot(
  "gtf_file does not exist" = file.exists(opts$gtf_file),
  "cytoband_file does not exist" = file.exists(opts$cytoband_file),
  "gene_order_file must end with `.txt.gz`" = stringr::str_ends(opts$gene_order_file, "\\.txt\\.gz")
)


# Read and prepare input files -----------------------------------------

gtf <- rtracklayer::import(opts$gtf_file, feature.type = "gene") |>
  as.data.frame() |>
  # rename to support joining
  dplyr::select(gene_id, chrom = seqnames, gene_start = start, gene_end = end) |>
  dplyr::mutate(chrom = glue::glue("chr{chrom}"))

cytoband_df <- readr::read_tsv(
  opts$cytoband_file,
  col_names = c("chrom", "chrom_arm_start", "chrom_arm_end", "band", "stain")
) |>
  # immediately remove chromosomes we don't want: contains `_`, or mito
  dplyr::filter(
    !stringr::str_detect(chrom, "_"),
    !stringr::str_starts(chrom, "chrM") # use M, not MT, to match both human and mouse
  ) |>
  # extract arm, if present, from band
  dplyr::mutate(arm = stringr::str_sub(band, 1, 1))

# Define chromosome factor orders ------------------------------------

# autosomes (numeric) chromosomes
autosome_chrom_df <- cytoband_df |>
  dplyr::mutate(chrom = stringr::str_remove(chrom, "chr")) |>
  dplyr::filter(!chrom %in% c("X", "Y")) |>
  dplyr::mutate(chrom = as.numeric(chrom)) |>
  dplyr::select(chrom, arm) |>
  dplyr::distinct() |>
  dplyr::arrange(chrom) |>
  dplyr::mutate(chrom_arm = glue::glue("{chrom}{arm}"))

# sex chromosomes
sex_chrom_df <- cytoband_df |>
  dplyr::mutate(chrom = stringr::str_remove(chrom, "chr")) |>
  dplyr::filter(chrom %in% c("X", "Y")) |>
  dplyr::select(chrom, arm) |>
  dplyr::distinct() |>
  dplyr::arrange(chrom) |>
  dplyr::mutate(chrom_arm = glue::glue("{chrom}{arm}"))

# combine autosome and sex dfs to derive final order
chrom_arm_levels <- glue::glue(
  "chr{c(autosome_chrom_df$chrom_arm, sex_chrom_df$chrom_arm)}"
)


# Prepare and export gene order file ------------------------------

gene_order_df <- cytoband_df |>
  # calculate arm boundaries
  dplyr::group_by(chrom, arm) |>
  dplyr::summarize(
    chrom_arm_start = min(chrom_arm_start),
    chrom_arm_end = max(chrom_arm_end),
    .groups = "drop"
  ) |>
  # combine gene coordinates with chromosome arm boundaries
  # use left_join to keep only cytoband_df chromosomes
  dplyr::left_join(
    gtf,
    by = "chrom",
    relationship = "many-to-many"
  ) |>
  # keep only genes actually on the given chromosome arm
  dplyr::filter(
    gene_start >= chrom_arm_start,
    gene_end <= chrom_arm_end
  ) |>
  # create chrom_arm column as identifier to use instead of chrom
  dplyr::mutate(chrom_arm = glue::glue("{chrom}{arm}")) |>
  # arrange and keep only relevant columns for infercnv
  dplyr::mutate(chrom_arm = factor(chrom_arm, levels = chrom_arm_levels)) |>
  dplyr::arrange(chrom_arm, gene_start) |>
  dplyr::select(gene_id, chrom_arm, gene_start, gene_end)

# Check that we got all the levels
stopifnot(
  "Missing or extra chromosomes in the gene order file" =
    setequal(
      chrom_arm_levels, levels(gene_order_df$chrom_arm)
    )
)

# export file without a header
readr::write_tsv(
  gene_order_df,
  opts$gene_order_file,
  col_names = FALSE
)
