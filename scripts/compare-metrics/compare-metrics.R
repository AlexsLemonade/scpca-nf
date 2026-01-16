#!/usr/bin/env Rscript

# get script location and activate renv
script_dir <- file.path(
  rprojroot::find_root(rprojroot::is_git_root),
  "scripts/compare-metrics"
)
renv::load(script_dir)

# parse options
library(optparse)

option_list <- list(
  make_option(
    c("-r", "--ref_s3"),
    type = "character",
    default = "s3://ccdl-scpca-results-prod-997241705947-us-east-1/results",
    help = "S3 URI for the bucket and prefix for the reference files"
  ),
  make_option(
    c("-c", "--comp_s3"),
    type = "character",
    default = "s3://ccdl-scpca-results-staging-997241705947-us-east-1/results",
    help = "S3 URI for the bucket and prefix for the comparison S3 files"
  ),
  make_option(
    c("-p", "--project_id"),
    type = "character",
    default = "",
    help = "Project ID(s) to filter the metrics files"
  ),
  make_option(
    c("-o", "--output_file"),
    type = "character",
    default = "",
    help = "Output file name for the rendered report. If not provided, a name will be generated with the current date and project ID(s)."
  ),
  make_option(
    c("--use_cache"), # temporary option for testing
    action = "store_true",
    default = FALSE,
    help = "Skip downloading files from S3 and use cached metrics data."
  )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# check parameters
stopifnot(
  "ref_s3 must be a valid S3 URI" = grepl("^s3://", opt$ref_s3),
  "comp_s3 must be a valid S3 URI" = grepl("^s3://", opt$comp_s3),
  "project_id must be provided" = nchar(opt$project_id) > 0
)

# split project ids (flexibly by commas, semicolons, and/or whitespace)
project_ids <- stringr::str_split_1(opt$project_id, "[,;\\s]+")

if (opt$output_file != "") {
  outfile <- opt$output_file
} else {
  timestamp <- format(Sys.time(), "%Y-%m-%d")
  outfile <- paste(
    timestamp,
    paste(project_ids, collapse = "-"),
    "metrics_comparison.html",
    sep = "_"
  )
}

# check that S3 paths are reachable
if (!(s3fs::s3_dir_exists(opt$ref_s3) && s3fs::s3_dir_exists(opt$comp_s3))) {
  stop(
    "Cannot access S3 location: ",
    opt$ref_s3,
    "and/or ",
    opt$comp_s3,
    "\n",
    "Please check that the paths exist and you have appropriate AWS credentials."
  )
}

rmarkdown::render(
  file.path(script_dir, "compare-metrics-template.rmd"),
  output_file = basename(outfile),
  output_dir = dirname(outfile),
  params = list(
    reference_s3 = opt$ref_s3,
    comparison_s3 = opt$comp_s3,
    project_id = project_ids,
    use_cache = opt$use_cache
  ),
  envir = new.env()
)
