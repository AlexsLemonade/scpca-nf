#!/usr/bin/env Rscript
#
# Run spell check and save results
# This script can be called by the pre-commit hook, in which case it should be
# given the changed files as arguments. Otherwise, it will check all R Markdown
# files in the repository.

arguments <- commandArgs(trailingOnly = TRUE)
file_pattern <- "\\.(Rmd|md|rmd)$"

# if there are arguments, check those files, otherwise check all markdown & rmd files
if (length(arguments) > 0) {
  precommit <- TRUE
  files <- arguments[grepl(file_pattern, arguments)]
} else {
  precommit <- FALSE
  # The only files we want to check are R Markdown and Markdown files
  files <- list.files(
    pattern = file_pattern,
    recursive = TRUE,
    full.names = TRUE
  )
}

# Find .git root directory
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# Read in dictionary
dict_file <- file.path(root_dir, "components", "dictionary.txt")
if (file.exists(dict_file)) {
  # Reading in the dictionary this way lets us put emojis in the dictionary file
  dictionary <- spelling::spell_check_files(dict_file)$word
} else {
  warning("Dictionary file not found")
  dictionary <- ""
}

# Run spell check
spelling_errors <- spelling::spell_check_files(files, ignore = dictionary) |>
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
