# This script takes the output folder from alevin-fry as input and 
# returns both the filtered and unfiltered counts matrices as a SingleCellExperiment stored in a .rds file 

# import libraries
library(magrittr)
library(optparse)
library(scpcaTools)

# set up arguments 
option_list <- list(
  make_option(
    opt_str = c("-r", "--run_id"),
    type = "character",
    help = "scpca_run_id"
  ), 
  make_option(
    opt_str = c("-c", "--which_counts"),
    type = "character", 
    help = "include counts for spliced cDNA only (spliced) or unspliced and spliced cDNA (unspliced)"
  ), 
  make_option(
    opt_str = c("-o", "-output_dir"),
    type = "character", 
    help = "output location for sce file"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

# check arguments are correct 
if(is.null(opt$run_id)){
  stop("Must input a scpca_run_id.")
}

# check for compatible sequencing unit types 
if(!(opt$which_counts %in% c("spliced", "unspliced"))){
  stop("Sequencing unit must be of type cell or nucleus")
}

# set up file directories 
# location where processed alevin-fry files live 
s3_processed_files <- 's3://nextflow-ccdl-scpca/scpca/processed'
sample_dir <- paste(file.path(s3_processed_files, opt$run_id), sep = "")

s3_output_files <- 's3://nextflow-ccdl-scpca/scpca/sce_processed'
sample_output_dir <- file.path(s3_output_files, opt$run_id)

unfiltered_sce_file <- file.path(sample_output_dir, "unfiltered_sce.rds")
filtered_sce_file <- file.path(sample_output_dir, "filtered_sce.rds")

if(!dir.exists(sample_output_dir)){
  dir.create(sample_output_dir, recursive = TRUE)
}


# get filtered and unfiltered sce 
# we may also want to eventually add in per cell QC metrics here 
unfiltered_sce <- import_quant_data(quant_dir = file.path(base_dir, sample), 
                                    tool = "alevin-fry", 
                                    which_counts = opt$which_counts, 
                                    usa_mode = TRUE)


filtered_sce <- import_quant_data(quant_dir = file.path(base_dir, sample), 
                                    tool = "alevin-fry", 
                                    which_counts = opt$which_counts, 
                                    usa_mode = TRUE, 
                                    filter = TRUE)

# write files 
readr::write_rds(unfiltered_sce, unfiltered_sce_file)
readr::write_rds(filtered_sce, filtered_sce_file)

