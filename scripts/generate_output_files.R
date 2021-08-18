# This script takes the output folder from alevin-fry as input and 
# returns both the filtered and unfiltered counts matrices as a SingleCellExperiment stored in a .rds file 

# import libraries
library(magrittr)
library(optparse)
library(scpcaTools)

# set up arguments 
option_list <- list(
  make_option(
    opt_str = c("-l", "--library_id"),
    type = "character",
    help = "scpca_library_id"
  ), 
  make_option(
    opt_str = c("-s", "--seq_unit"),
    type = "character", 
    help = "include counts for spliced cDNA only (spliced) or unspliced and spliced cDNA (unspliced)"
  ), 
  make_option(
    opt_str = c("-o", "--output_dir"),
    type = "character", 
    help = "output location for sce file"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

# check arguments are correct 
if(is.null(opt$library_id)){
  stop("Must input a scpca_library_id.")
}

# check for compatible sequencing unit types 
if(!(opt$seq_unit %in% c("cell", "nucleus"))){
  stop("Sequencing unit must be of type cell or nucleus")
}

# set up filenames
unfiltered_sce_file <- file.path(opt$output_dir, "unfiltered_sce.rds")
filtered_sce_file <- file.path(opt$output_dir, "filtered_sce.rds")

which_counts <- dplyr::case_when(opt$seq_unit == "cell" ~ "spliced",
                                 opt$seq_unit == "nucleus" ~ "unspliced")

# get filtered and unfiltered sce 
unfiltered_sce <- import_quant_data(quant_dir = opt$output_dir, 
                                    tool = "alevin-fry", 
                                    which_counts = opt$seq, 
                                    usa_mode = TRUE)


filtered_sce <- import_quant_data(quant_dir = opt$output_dir, 
                                    tool = "alevin-fry", 
                                    which_counts = opt$which_counts, 
                                    usa_mode = TRUE, 
                                    filter = TRUE)

# write files 
readr::write_rds(unfiltered_sce, unfiltered_sce_file)
readr::write_rds(filtered_sce, filtered_sce_file)

