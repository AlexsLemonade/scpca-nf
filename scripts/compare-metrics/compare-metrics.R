#!/usr/bin/env Rscript

# get script location and activate renv
script_dir <- here::here("scripts/compare-metrics")
renv::load(script_dir)

# parse options
library(optparse)
