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
      Must contain the `organ` column, `--cell_type_column` and `--gene_id_column`. "
  ),
  make_option(
    opt_str = c("--cell_type_column"),
    type = "character",
    default = "cell type",
    help = "Column name in `marker_gene_file` containing cell types."
  ),
  make_option(
    opt_str = c("--gene_id_column"),
    type = "character",
    default = "official gene symbol",
    help = "Column name in `marker_gene_file` containing the gene symbols."
  ),
  make_option(
    opt_str = c("--gtf_file"),
    type = "character",
    help = "reference gtf to use"
  ),
  make_option(
    opt_str = c("--ref_mtx_file"),
    type = "character",
    help = "path to save reference binary matrix file"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

# Set up -----------------------------------------------------------------------

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

# Create marker gene matrix ----------------------------------------------------

# grab marker genes for specified organs from marker gene reference
organ_gene_df <- marker_gene_df |>
  dplyr::filter(organ %in% organs) |>
  dplyr::select({{cell_type_column}}, {{gene_id_column}})

# create a binary matrix of genes by cell type markers 
binary_mtx <- organ_gene_df |>
  tidyr::pivot_wider(id_cols = gene_id_column,
                     names_from = cell_type_column,
                     values_from = cell_type_column,
                     values_fn = length,
                     values_fill = 0) |> 
  tibble::column_to_rownames(gene_id_column) |>
  # add a column with no marker genes 
  # cell assign will assign cells to "other" when no other cell types are appropriate 
  dplyr::mutate(other = 0)

# replace length with 1
binary_mtx[binary_mtx > 1] <- 1

# Add ensembl ids --------------------------------------------------------------

# combine with ensembl gene id 
# read in gtf file (genes only for speed)
gtf <- rtracklayer::import(opt$gtf_file, feature.type = "gene")

# create a data frame with ensembl id and gene symbol 
gene_id_map <- gtf |>
  as.data.frame() |>
  dplyr::select(
    "ensembl_id" = "gene_id",
    "gene_symbol" = "gene_name"
  ) |>
  dplyr::filter(gene_symbol != "NA") |>
  dplyr::distinct() |>
  ## in case there are any duplicate gene_ids (there shouldn't be!)
  dplyr::group_by(ensembl_id) |>
  dplyr::summarize(gene_symbol = paste(gene_symbol, collapse = ";"))

# replace gene symbols with ensembl id
binary_mtx <- binary_mtx |>
  tibble::rownames_to_column("gene_symbol") |>
  dplyr::left_join(gene_id_map, by = "gene_symbol") |>
  # drop any rows where there is no ensembl id
  tidyr::drop_na(ensembl_id) |>
  dplyr::select(ensembl_id, everything(), -gene_symbol)

readr::write_tsv(binary_mtx, opt$ref_mtx_file)
