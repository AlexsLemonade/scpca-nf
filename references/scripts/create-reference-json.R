#!/usr/bin/env Rscript

# Script for creating json file to hold path to reference files for organisms that can be
# used with scpca-nf. The output of this script will create a json file where each key corresponds
# to an organism and contains a dictionary of reference paths. The paths included here are specific
# to the organization used for storing references in `s3://scpca-references`.
# To create the json file, use a TSV file that contains three columns, `organism`, `assembly`, and `version`

library(optparse)
project_root <- rprojroot::find_root(rprojroot::has_dir(".git"))

option_list <- list(
  make_option(
    opt_str = c("--ref_metadata"),
    type = "character",
    default = file.path(project_root, "references", "ref-metadata.tsv"),
    help = ""
  ),
  make_option(
    opt_str = c("--ref_json"),
    type = "character",
    default = file.path(project_root, "references", "scpca-refs.json")
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

# Function for creating json entry ---------------------------------------------

# function to create individual json entries containing individual reference paths
create_ref_entry <- function(
    organism,
    assembly,
    version,
    reference_name) {
  # create base reference directory
  ref_dir <- file.path(
    tolower(organism),
    glue::glue("ensembl-{version}")
  )
  fasta_dir <- file.path(ref_dir, "fasta")
  annotation_dir <- file.path(ref_dir, "annotation")

  # create a single json entry containing all necessary file paths
  json_entry <- list(
    ref_dir = ref_dir,
    ref_fasta = file.path(
      fasta_dir,
      glue::glue("{organism}.{assembly}.dna.primary_assembly.fa.gz")
    ),
    ref_fasta_index = file.path(
      fasta_dir,
      glue::glue("{organism}.{assembly}.dna.primary_assembly.fa.fai")
    ),
    ref_gtf = file.path(
      annotation_dir,
      glue::glue("{reference_name}.gtf.gz")
    ),
    mito_file = file.path(
      annotation_dir,
      glue::glue("{reference_name}.mitogenes.txt")
    ),
    t2g_3col_path = file.path(
      annotation_dir,
      glue::glue("{reference_name}.spliced_intron.tx2gene_3col.tsv")
    ),
    t2g_bulk_path = file.path(
      annotation_dir,
      glue::glue("{reference_name}.spliced_cdna.tx2gene.tsv")
    ),
    splici_index = file.path(
      ref_dir, "salmon_index",
      glue::glue("{reference_name}.spliced_intron.txome")
    ),
    salmon_bulk_index = file.path(
      ref_dir, "salmon_index",
      glue::glue("{reference_name}.spliced_cdna.txome")
    ),
    cellranger_index = file.path(
      ref_dir, "cellranger_index",
      glue::glue("{reference_name}_cellranger_full")
    ),
    star_index = file.path(
      ref_dir, "star_index",
      glue::glue("{reference_name}.star_idx")
    )
  )

  return(json_entry)
}

# Compile json file ------------------------------------------------------------

# check that input metadata exists and read in
if (!file.exists(opt$ref_metadata)) {
  stop("ref_metadata file does not exist.")
}

ref_metadata <- readr::read_tsv(opt$ref_metadata, col_types = "c") |>
  dplyr::mutate(reference_name = glue::glue("{organism}.{assembly}.{version}"))

# get entries for all organisms provided
all_entries <- ref_metadata |>
  purrr::pmap(create_ref_entry) |>
  purrr::set_names(ref_metadata$reference_name)

# Write to JSON
jsonlite::write_json(
  all_entries,
  path = opt$ref_json,
  pretty = TRUE,
  auto_unbox = TRUE
)
