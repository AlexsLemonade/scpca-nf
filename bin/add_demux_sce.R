#!/usr/bin/env Rscript

# import libraries
suppressPackageStartupMessages({
  library(optparse)
  library(SingleCellExperiment)
})
# set up arguments
option_list <- list(
  make_option(
    opt_str = c("-s", "--sce_file"),
    type = "character",
    help = "path to rds file with an sce object"
  ),
  make_option(
    opt_str = c("-l", "--library_id"),
    type = "character",
    help = "library id for the sce object; required for filtering cellhash pool file",
    default = NULL
  ),
  make_option(
    opt_str = c("-o", "--output_sce_file"),
    type = "character",
    help = "path to output with multiplex data added. Must end in .rds"
  ),
  make_option(
    opt_str = c("-v", "--vireo_dir"),
    type = "character",
    help = "path to vireo output directory",
    default = NULL
  ),
  make_option(
    opt_str = c("--cellhash_pool_file"),
    type = "character",
    help = "path to table of cellhash barcodes and sample ids",
    default = NULL
  ),
  make_option(
    opt_str = c("-H", "--hash_demux"),
    action = "store_true",
    help = "add HashedDrops demultiplex results to sce"
  ),
  make_option(
    opt_str = c("-S", "--seurat_demux"),
    action = "store_true",
    help = "add Seurat demultiplex results to sce"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))


# check that unfiltered file file exists
if(!file.exists(opt$sce_file)){
  stop("Can't find input SCE file")
}

# check that output file name ends in .rds
if(!(stringr::str_ends(opt$output_sce_file, ".rds"))){
  stop("output file name must end in .rds")
}

# check for donor_ids file in vireo_dir
if(!is.null(opt$vireo_dir)){
  vireo_file <- file.path(opt$vireo_dir, "donor_ids.tsv")
  if(!file.exists(vireo_file)){
    stop("Can't find donor_ids.tsv file in vireo directory")
  }
}

if(!is.null(opt$cellhash_pool_file)){
  if(!file.exists(opt$cellhash_pool_file)){
    stop("Can't find cellhash_pool_file")
  }
  if(is.null(opt$library_id)){
    stop("Must specify library_id with cellhash_pool_file")
  }
  cellhash_df <- readr::read_tsv(opt$cellhash_pool_file) |>
    dplyr::filter(scpca_library_id == opt$library_id) |>
    dplyr::select(sample_id = scpca_sample_id, barcode_id)
}

# read in sce rds file
sce <- readRDS(opt$sce_file)

# check for cellhash altExp if we will use it
if(!is.null(cellhash_df) || opt$hash_demux || opts$seurat_demux){
  if(!"cellhash" %in% altExpNames(sce)){
    stop("Can't process cellhash demulitplexing without a 'cellhash' altExp")
  }
}

# add cellhash sample data to SCE
if(!is.null(cellhash_df)){
  sce <- scpcaTools::add_cellhash_ids(sce, cellhash_df, remove_unlabeled = TRUE)
}

# add HashedDrops results
if(opt$hash_demux){
  sce <- scpcaTools::add_demux_hashedDrops(sce)
}

# add Seurat results
if(opt$seurat_demux){
  sce <- scpcaTools::add_demux_seurat(sce)
}

# add vireo results
if(!is.null(opt$vireo_dir)){
  vireo_table <- readr::read_tsv(vireo_file)
  sce <- scpcaTools::add_demux_vireo(sce, vireo_table)
}



# write filtered sce to output
readr::write_rds(sce, opt$output_sce_file, compress = "gz")
