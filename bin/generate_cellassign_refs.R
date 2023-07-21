#!/usr/bin/env Rscript

# this script is used to create binary gene x cell marker genes matrices 
# for use as references with CellAssign

library(optparse)

option_list <- list(
  make_option(
    opt_str = c("--organs"),
    type = "character",
    help = "Comma separated list of organs to include in gene x cell matrix"
  ),
  make_option(
    opt_str = c("--marker_gene_file"),
    type = "character",
    help = "File containing a list of marker genes for all cell types and all organs. 
      Must contain the `organ`, `cell_type_column` and `gene_id_column` columns. "
  ),
  make_option(
    opt_str = c("--cell_type_column"),
    type = "character",
    default = "cell type"
  ),
  make_option(
    opt_str = c("--gene_id_column"),
    type = "character",
    default = "official gene symbol"
  ),
  make_option(
    opt_str = c("--gtf"),
    type = "character",
    help = "reference gtf to use"
  ),
  make_option(
    opt_str = c("--ref_file"),
    type = "character",
    help = "path to save reference file"
  )
)


if(!file.exists(opt$marker_gene_file)){
  stop("Provided `marker_gene_file` does not exist.")
}

marker_gene_df <- readr::read_tsv(opt$marker_gene_file)

organs <- stringr::str_split(opt$organs, ",") |>
  unlist() |>
  stringr::str_trim()

if(!all(organs %in% marker_gene_df$organ)){
  stop("`organs` must be present in `organ` column of `marker_gene_file`.")
}

cell_type_column <- opt$cell_type_column
gene_id_column <- opt$gene_id_column

if(!(cell_type_column %in% colnames(marker_gene_df))){
  stop("Specified `cell_type_column` does not exist as a column in `marker_gene_file`")
}

if(!(gene_id_column %in% colnames(marker_gene_df))){
  stop("Specified `gene_id_column` does not exist as a column in `marker_gene_file`")
}

organ_gene_df <- marker_gene_df |>
  dplyr::filter(organ %in% organs) |>
  dplyr::select(cell_type_column, gene_id_column)

# combine with ensembl gene id 

# create a binary matrix of genes by cell type markers 
binary_mtx <- organ_gene_df |>
  tidyr::pivot_wider(id_cols = ensembl_id,
                     names_from = celltype_column,
                     values_from = celltype_column,
                     values_fn = length,
                     values_fill = 0) |> 
  tibble::column_to_rownames(ensembl_id) |>
  # add a column with no marker genes 
  # cell assign will assign cells to "other" when no other cell types are appropriate 
  dplyr::mutate(other = 0)

# replace length with 1
binary_mtx[binary_mtx > 1] <- 1
