#!/usr/bin/env Rscript

# import libraries
library(optparse)

# Set up optparse options
option_list <- list(
  make_option(
    opt_str = c("--reference_files"),
    type = "character",
    default = NULL,
    help = "comma-separated list of reference file names"
  ),
  make_option(
    opt_str = c("-o", "--output_file"),
    type = "character",
    default = NULL,
    help = "path to output tsv file"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

ref_files <- stringr::str_split_1(opt$reference_files, ",") |>
  basename()

ref_info <- data.frame(filename = ref_files) |>
  dplyr::mutate(
    ref_string = tools::file_path_sans_ext(filename)
  ) |>
  tidyr::separate_wider_delim(
    ref_string,
    names = c("reference_name", "source", "version"),
    delim = "_",
    too_many = "drop"
  ) |>
  dplyr::mutate(
    # replace `-` with `.` for non-date versions
    version = ifelse(
      grepl("^20[0-9]{2}-[0-9]{2}-[0-9]{2}$", version),
      version,
      stringr::str_replace_all(version, "-", ".")
    )
  )
