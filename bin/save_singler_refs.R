#!/usr/bin/env Rscript

# this script is used to download celldex references and save them as RDS files

library(optparse)
library(celldex)

option_list <- list(
  make_option(
    opt_str = c("--ref_name"),
    type = "character",
    help = "name associated with celldex reference, e.g., HumanPrimaryCellAtlasData"
  ),
  make_option(
    opt_str = c("--ref_file"),
    type = "character",
    help = "path to save reference file"
  )
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

# check that provided ref name is in celldex package
if(! opt$ref_name %in% ls("package:celldex")){
  stop(glue::glue("Provided `ref_name` `{opt$ref_name}` does not match a celldex dataset."))
} 

# get a reference library from celldex:
ref <- do.call(`::`, args = list("celldex", opt$ref_name))(ensembl = TRUE)

# export ref data 
readr::write_rds(ref, opt$ref_file, compress = "gz")
