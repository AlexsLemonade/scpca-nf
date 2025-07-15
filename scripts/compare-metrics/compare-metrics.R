#!/usr/bin/env Rscript

# get script location and activate renv
script_dir <- file.path(rprojroot::find_root(rprojroot::is_git_root), "scripts/compare-metrics")
renv::load(script_dir)

# parse options
library(optparse)

option_list <- list(
  make_option(
    c("-r", "--ref_s3"),
    type = "character",
    default = "s3://nextflow-ccdl-results/scpca-prod/results",
    help = "S3 URI for the bucket and prefix for the reference files"
  ),
  make_option(
    c("-c", "--comp_s3"),
    type = "character",
    default = "s3://nextflow-ccdl-results/metrics-test", # temporary default for testing
    # default = "s3://nextflow-ccdl-results/scpca-staging/results", # normal default
    help = "S3 URI for the bucket and prefix for the comparison S3 files"
  ),
  make_option(c("-p", "--project_id"),
    type = "character",
    default = "all",
    help = "Project ID(s) to filter the metrics files. Default is all projects"
  ),
  make_option(
    c("-o", "--output_file"),
    type = "character",
    default = "compare-metrics.html",
    help = "Output file name for the rendered report"
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
  "comp_s3 must be a valid S3 URI" = grepl("^s3://", opt$comp_s3)
)

# split project ids (flexibly by commas, semicolons, and/or whitespace)
project_ids <- stringr::str_split_1(opt$project_id, "[,;\\s]+")

rmarkdown::render(
  file.path(script_dir, "compare-metrics-template.rmd"),
  output_file = basename(opt$output_file),
  output_dir = dirname(opt$output_file),
  params = list(
    reference_s3 = opt$ref_s3,
    comparison_s3 = opt$comp_s3,
    project_id = project_ids,
    use_cache = opt$use_cache
  ),
  envir = new.env()
)
