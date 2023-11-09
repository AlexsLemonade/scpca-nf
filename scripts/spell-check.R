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
dictionary <- readLines(dict_file, encoding = "UTF-8")

# Add emoji to dictionary
dictionary_plus <- c(dictionary, spelling::spell_check_text("⚠️")$word)

# The only files we want to check are R Markdown and Markdown files
files <- list.files(pattern = "\\.(Rmd|md|rmd)$", recursive = TRUE, full.names = TRUE)

# Run spell check
spelling_errors <- spelling::spell_check_files(files, ignore = dictionary_plus)
spelling_tibble <- spelling_errors |>
  data.frame() |>
  tidyr::unnest(cols = found) |>
  tidyr::separate(found, into = c("file", "lines"), sep = ":") |>
  dplyr::distinct()


if (precommit) {
  if (nrow(spelling_errors) > 0) {
    cat("The following spelling errors were found:\n")
    print(spelling_errors)
  }

  # Update dictionary for future use
  updated_dict <- union(dictionary, spelling_errors$word)
  # remove empty strings
  updated_dict <- updated_dict[updated_dict != ""]
  writeLines(sort(updated_dict), dict_file)
} else {
  # Print out how many spell check errors
  write(nrow(spelling_tibble), stdout())
  # Save spell errors to file temporarily
  readr::write_tsv(spelling_tibble, "spell_check_errors.tsv")
}
