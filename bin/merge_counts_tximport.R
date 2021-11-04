#!/usr/bin/env Rscript

# This script grabs the quant.sf output files from Salmon for all samples within a project
# The files are then imported use tximport to generate a gene x count matrix
# Output includes a counts matrix as a tab separated tsv file and a .Rds object containing the tximport object


# load needed packages
library(magrittr)
library(tximport)
library(optparse)

# Set up optparse options
option_list <- list(
  make_option(
    opt_str = c("-p", "--project_id"),
    type = "character",
    help = "scpca project ID",
  ),
  make_option(
    opt_str = c("-s", "--salmon_dir"),
    type = "character",
    help = "Path to salmon output directories"
  ),
  make_option(
    opt_str = c("-o", "--output_file"),
    type = "character",
    help = "Path to output RDS file, must end in .rds"
  ),
  make_option(
    opt_str = c("-t", "--tx2gene"),
    type = "character",
    help = "Path to tx2gene file created from index used for mapping."
  )
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

# read in tx2gene
tx2gene <- readr::read_tsv(opt$tx2gene, col_names = c("transcript_id", "gene_id"))

# list of paths to salmon files 
library_ids <- readLines(opt$salmon_dir)
salmon_files <- file.path(library_ids, "quant.sf")
names(salmon_files) <- library_ids

# import using tximport
txi_salmon <- tximport(salmon_files, type = "salmon", tx2gene = tx2gene)

# write counts matrix to txt file
txi_salmon$counts %>%
  as.data.frame %>%
  tibble::rownames_to_column("gene_id") %>%
  readr::write_tsv(file = opt$output_file)
