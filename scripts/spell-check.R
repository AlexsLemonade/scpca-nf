#!/usr/bin/env Rscript
#
# Run spell check and save results
# Adapted from: https://github.com/AlexsLemonade/refinebio-examples/blob/33cdeff66d57f9fe8ee4fcb5156aea4ac2dce07f/scripts/spell-check.R

# Find .git root directory
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# Read in dictionary
dict_file <- file.path(root_dir, "components", "dictionary.txt")
dictionary <- readLines(dict_file)

# Add emoji to dictionary
dictionary <- c(dictionary, spelling::spell_check_text("⚠️")$word)

# The only files we want to check are R Markdown and Markdown files
files <- list.files(pattern = "\\.(Rmd|md|rmd)$", recursive = TRUE, full.names = TRUE)

# Run spell check
spelling_errors <- spelling::spell_check_files(files, ignore = dictionary) |>
  data.frame() |>
  tidyr::unnest(cols = found) |>
  tidyr::separate(found, into = c("file", "lines"), sep = ":")

# Print out how many spell check errors
write(nrow(spelling_errors), stdout())

# Save spell errors to file temporarily
readr::write_tsv(spelling_errors, "spell_check_errors.tsv")

# Update dictionary for future use
dictionary <- sort(unique(c(dictionary, spelling_errors$word)))
writeLines(dictionary, dict_file)
