#!/usr/bin/env Rscript
#
# Run spell check and save results
# Adapted from: https://github.com/AlexsLemonade/refinebio-examples/blob/33cdeff66d57f9fe8ee4fcb5156aea4ac2dce07f/scripts/spell-check.R

args <- commandArgs(trailingOnly = TRUE)
precommit <- "--precommit" %in% args

# Find .git root directory
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# Read in dictionary
dict_file <- file.path(root_dir, "components", "dictionary.txt")
dictionary <- readLines(dict_file)

# Add emoji to dictionary
dictionary_plus <- c(dictionary, spelling::spell_check_text("⚠️")$word)

# The only files we want to check are R Markdown and Markdown files
files <- list.files(pattern = "\\.(Rmd|md|rmd)$", recursive = TRUE, full.names = TRUE)
# remove any files in the work directory
files <- files[!grepl("^./work/", files)]

# Run spell check
spelling_errors <- spelling::spell_check_files(files, ignore = dictionary_plus) |>
  data.frame() |>
  tidyr::unnest(cols = found) |>
  tidyr::separate(found, into = c("file", "lines"), sep = ":")

if (precommit) {
  if (nrow(spelling_errors) > 0) {
    cat("The following spelling errors were found:\n")
    print(data.frame(spelling_errors))
  }

  # Update dictionary for future use
  updated_dict <- union(dictionary, spelling_errors$word)
  # remove empty strings
  updated_dict <- updated_dict[updated_dict != ""]
  # case insensitive sort
  updated_dict <- updated_dict[order(tolower(updated_dict))]
  writeLines(updated_dict, dict_file)
} else {
  # Print out how many spell check errors
  write(nrow(spelling_errors), stdout())
  # Save spell errors to file temporarily
  readr::write_tsv(spelling_errors, "spell_check_errors.tsv")
}
