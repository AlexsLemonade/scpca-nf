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
  ),
  make_option(
    opt_str = c("-t", "--technology"),
    type = "character",
    help = "sequencing technology string to store in metadata"
  ),
  make_option(
    opt_str = c("--library_id"),
    type = "character",
    help = "library id"
  ),
  make_option(
    opt_str = c("--sample_id"),
    type = "character",
    help = "sample id(s). If more than one, separated by commas and/or semicolons."
  ),
  make_option(
    opt_str = c("--sample_metadata_file"),
    type = "character",
    help = "path to tsv file containing sample metadata"
  ),
  make_option(
    opt_str = c("--spliced_only"),
    action = "store_true",
    default = FALSE,
    help = "include only the spliced counts as the main counts assay in the returned SCE object"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

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

# check that sample metadata file exists
if(!file.exists(opt$sample_metadata_file)){
  stop("sample metadata file not found.")
}

# read in mitochondrial gene list
mito_genes <- unique(scan(opt$mito_file, what = "character"))

# read in gtf file (genes only for speed)
gtf <- rtracklayer::import(opt$gtf_file, feature.type = "gene")

# parse sample id list
sample_ids <- unlist(stringr::str_split(opt$sample_id, ",|;")) |> sort()

# read in sample metadata and filter to sample ids
sample_metadata_df <- readr::read_tsv(opt$sample_metadata_file) |>
  dplyr::filter(scpca_sample_id %in% sample_ids)

# set include unspliced for non feature data
include_unspliced <- !opt$spliced_only

# get unfiltered sce
unfiltered_sce <- read_alevin(quant_dir = opt$alevin_dir,
                              include_unspliced = include_unspliced,
                              fry_mode = TRUE,
                              tech_version = opt$technology,
                              library_id = opt$library_id,
                              sample_id = sample_ids)


# read and merge feature counts if present
if (opt$feature_dir != ""){
  feature_sce <- read_alevin(quant_dir = opt$feature_dir,
                             include_unspliced = FALSE,
                             fry_mode = TRUE,
                             feature_data = TRUE,
                             library_id = opt$library_id,
                             sample_id = sample_ids)

  unfiltered_sce <- merge_altexp(unfiltered_sce, feature_sce, opt$feature_name)
  # add alt experiment features stats
  altExp(unfiltered_sce, opt$feature_name) <- scuttle::addPerFeatureQCMetrics(altExp(unfiltered_sce, opt$feature_name))
}

# add per cell and per gene statistics to colData and rowData
unfiltered_sce <- unfiltered_sce |>
  add_cell_mito_qc(mito = mito_genes) |>
 # add gene symbols to rowData
  add_gene_symbols(gene_info = gtf) |>
  scuttle::addPerFeatureQCMetrics()

# get a list of sample metadata to add
# for each sample in the library, create an individual list of sample metadata
sample_metadata_list <- purrr::map(sample_ids,
                                   \(sample){
                                     single_sample_df <- sample_metadata_df |>
                                       dplyr::filter(scpca_sample_id %in% sample)
                                     list(
                                       age = single_sample_df$age,
                                       sex = single_sample_df$sex,
                                       diagnosis = single_sample_df$sex,
                                       subdiagnosis = single_sample_df$subdiagnosis,
                                       tissue_location = single_sample_df$tissue_location,
                                       disease_timing = single_sample_df$disease_timing,
                                       organism = single_sample_df$organism,
                                       development_stage_ontology_term_id = single_sample_df$development_stage_ontology_term_id,
                                       sex_ontology_term_id = single_sample_df$sex_ontology_term_id,
                                       organism_ontology_id = single_sample_df$organism_ontology_id,
                                       self_reported_ethnicity_ontology_term_id = single_sample_df$self_reported_ethnicity_ontology_term_id,
                                       disease_ontology_term_id = single_sample_df$disease_ontology_term_id,
                                       tissue_ontology_term_id = single_sample_df$tissue_ontology_term_id
                                     )
                                   }) |>
  purrr::set_names(sample_ids)

# add list of sample metadata
# access the list of sample metadata for each sample
# with metadata(unfiltered_sce)$sample_metadata[[sample_id]]
metadata(unfiltered_sce)$sample_metadata <- sample_metadata_list

# write to rds
readr::write_rds(unfiltered_sce, opt$unfiltered_file, compress = "gz")

