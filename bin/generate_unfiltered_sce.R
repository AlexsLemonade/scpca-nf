#!/usr/bin/env Rscript

# This script takes the output folder from alevin-fry as input and
# returns the unfiltered counts matrices as a SingleCellExperiment stored in a .rds file

# import libraries
suppressPackageStartupMessages({
  library(optparse)
  library(SingleCellExperiment)
  library(scpcaTools)
})
# set up arguments
option_list <- list(
  make_option(
    opt_str = c("-s", "--seq_unit"),
    type = "character",
    help = "`cell` or `nucleus`, which will determine whether to include counts for spliced cDNA only (cell) or unspliced and spliced cDNA (nucleus)"
  ),
  make_option(
    opt_str = c("-a", "--alevin_dir"),
    type = "character",
    help = "directory with alevin output files for RNA-seq quantification"
  ),
  make_option(
    opt_str = c("-f", "--feature_dir"),
    type = "character",
    default = "",
    help = "directory with alevin output files for feature quantification"
  ),
  make_option(
    opt_str = c("-n", "--feature_name"),
    type = "character",
    default = "ALT",
    help = "Feature type"
  ),
  make_option(
    opt_str = c("-u", "--unfiltered_file"),
    type = "character",
    help = "path to output unfiltered rds file. Must end in .rds"
  ),
  make_option(
    opt_str = c("-m", "--mito_file"),
    type = "character",
    help = "path to list of mitochondrial genes"
  ),
  make_option(
    opt_str = c("-g", "--gtf_file"),
    type = "character",
    help = "path to gtf file with gene annotations"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

# check for compatible sequencing unit types
if(!(opt$seq_unit %in% c("cell", "nucleus"))){
  stop("Sequencing unit must be of type cell or nucleus")
}

# check that output file name ends in .rds
if(!(stringr::str_ends(opt$unfiltered_file, ".rds"))){
  stop("unfiltered file name must end in .rds")
}

# check that mitochondrial gene list exists
if(!file.exists(opt$mito_file)){
  stop("Mitochondrial gene list file not found.")
}

# check that gtf file exists
if(!file.exists(opt$gtf_file)){
  stop("gtf file not found.")
}

# read in mitochondrial gene list
mito_genes <- unique(scan(opt$mito_file, what = "character"))

# read in gtf file (genes only for speed)
gtf <- rtracklayer::import(opt$gtf_file, feature_type = "gene")

# convert seq_unit to spliced or unspliced to determine which types of transcripts to include in final counts matrix
which_counts <- dplyr::case_when(opt$seq_unit == "cell" ~ "spliced",
                                 opt$seq_unit == "nucleus" ~ "unspliced")

# get unfiltered sce
unfiltered_sce <- read_alevin(quant_dir = opt$alevin_dir,
                              which_counts = which_counts,
                              usa_mode = TRUE)


# read and merge feature counts if present
if (opt$feature_dir != ""){
  feature_sce <- read_alevin(quant_dir = opt$feature_dir,
                             mtx_format = TRUE) 
   
  unfiltered_sce <- merge_altexp(unfiltered_sce, feature_sce, opt$feature_name) 
  # add alt experiment features stats
  altExp(unfiltered_sce, opt$feature_name) <- scuttle::addPerFeatureQCMetrics(altExp(unfiltered_sce, opt$feature_name))
}

# add per cell and per gene statistics to colData and rowData
unfiltered_sce <- unfiltered_sce |>
  add_cell_mito_qc(mito = mito_genes) |>
  add_gene_symbols(gene_info = gtf) |>
  scuttle::addPerFeatureQCMetrics()

# write to rds
readr::write_rds(unfiltered_sce, opt$unfiltered_file, compress = "gz")

