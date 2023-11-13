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
    opt_str = c("--ref_file_prefix"),
    type = "character",
    help = "prefix to use for naming saved reference file.
      All files will be saved with the following name: <ref_file_prefix>_<celldex_version>.rds"
  )
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

# get celldex package version
celldex_version <- packageVersion("celldex") |>
  as.character() |>
  # replace . in version string to avoid . in filename
  stringr::str_replace_all("\\.", "-")

# construct output ref file name
ref_output_file <- glue::glue("{opt$ref_file_prefix}_{celldex_version}.rds")

# check that provided ref name is in celldex package
if (!opt$ref_name %in% ls(getNamespace("celldex"))) {
  stop(glue::glue("Provided `ref_name` `{opt$ref_name}` does not match a celldex dataset."))
}

# get a reference library from celldex:
ref <- do.call(`::`, args = list("celldex", opt$ref_name))(ensembl = TRUE)

# export ref data
readr::write_rds(ref, ref_output_file, compress = "gz")
