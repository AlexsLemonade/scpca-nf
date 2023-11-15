#!/usr/bin/env Rscript

# A small script to read in a comma separated list of celltype reference files
# and output a table with the reference name, source, and version

# File names must be of the format:
#  <reference_name>_<source>_<version><optional_extra.ext>

# Usage:
# Rscript make_celltype_ref_table.R <list,of,file,paths> <output_table.tsv>



# not using optparse as to avoid dependencies
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Incorrect number of arguments")
}

reference_files <- args[1]
output_file <- args[2]

if (!stringr::str_ends(output_file, ".tsv")) {
  stop("Output file must end in .tsv")
}

ref_files <- stringr::str_split_1(reference_files, ",") |>
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

readr::write_tsv(ref_info, output_file)
