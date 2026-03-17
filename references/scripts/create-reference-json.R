#!/usr/bin/env Rscript

# Script for creating json file to hold path to reference files for organisms that can be
# used with scpca-nf. The output of this script will create a json file where each key corresponds
# to an organism and contains a dictionary of reference paths. The paths included here are specific
# to the organization used for storing references in `s3://scpca-nf-references`.
# To create the json file, use a TSV file that contains the following columns:
# `organism`, `assembly`, `version`, `include_salmon`, `include_cellranger`, `include_star`

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
  include_salmon,
  include_cellranger,
  include_star,
  include_infercnv,
  include_flex,
  include_visium,
  reference_name
) {
  # create base reference directory
  ref_dir <- file.path(
    tolower(organism),
    glue::glue("ensembl-{version}")
  )
  fasta_dir <- file.path(ref_dir, "fasta")
  annotation_dir <- file.path(ref_dir, "annotation")
  flex_probe_dir <- file.path(ref_dir, "flex-probe-refs")
  visium_probe_dir <- file.path(ref_dir, "visium-probe-refs")
  # cell ranger uses these values in their file names rather than the proper organism name
  cellranger_organism_map <- c(
    "Homo_sapiens" = "Human",
    "Mus_musculus" = "Mouse"
  )
  cellranger_organism <- cellranger_organism_map[[organism]]

  # map for visium probe files to ensure we only record the relevant ones
  visium_probe_map <- tibble::tribble(
    ~technology, ~organism, ~internal_assembly, ~ensembl_version, ~visium_assembly,
    "visium1_v1", "Human", "GRCh38", "98", "GRCh38",
    "visium2_v2", "Human", "GRCh38", "98", "GRCh38",
    "visium2_v2.1", "Human", "GRCh38", "110", "GRCh38",
    "visium1_v1", "Mouse", "GRCm38", "98", "mm10",
    "visium2_v2", "Mouse", "GRCm38", "98", "mm10",
    "visium2_v2.1", "Mouse", "GRCm39", "110", "GRCm39"
  )

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
    # fill in optional entries with empty strings
    t2g_3col_path = "",
    t2g_bulk_path = "",
    splici_index = "",
    salmon_bulk_index = "",
    cellranger_index = "",
    star_index = "",
    infercnv_gene_order = "",
    cytoband = "",
    flex_probe_files = "",
    visium_probe_files = ""
  )

  # fill in values related to salmon/alevin-fry index
  if (include_salmon) {
    json_entry$t2g_3col_path <- file.path(
      annotation_dir,
      glue::glue("{reference_name}.spliced_intron.tx2gene_3col.tsv")
    )
    json_entry$t2g_bulk_path <- file.path(
      annotation_dir,
      glue::glue("{reference_name}.spliced_cdna.tx2gene.tsv")
    )
    json_entry$splici_index <- file.path(
      ref_dir,
      "salmon_index",
      glue::glue("{reference_name}.spliced_intron.txome")
    )
    json_entry$salmon_bulk_index <- file.path(
      ref_dir,
      "salmon_index",
      glue::glue("{reference_name}.spliced_cdna.txome")
    )
  }

  # fill in values related to cellranger index
  if (include_cellranger) {
    json_entry$cellranger_index <- file.path(
      ref_dir,
      "cellranger_index",
      glue::glue("{reference_name}_cellranger_full")
    )
  }

  # fill in values related to star index
  if (include_star) {
    json_entry$star_index <- file.path(
      ref_dir,
      "star_index",
      glue::glue("{reference_name}.star_idx")
    )
  }

  # fill in values related to infercnv gene order file
  if (include_infercnv) {
    json_entry$infercnv_gene_order <- file.path(
      ref_dir,
      "infercnv",
      glue::glue("{reference_name}_gene_order_arms.txt.gz")
    )
    json_entry$cytoband <- file.path(
      annotation_dir,
      glue::glue("{reference_name}_cytoband.txt.gz")
    )
  }

  # add directory for flex probes
  if (include_flex) {
    json_entry$flex_probe_files <- list(
      "10xflex_v1.1_single" = file.path(
        flex_probe_dir,
        glue::glue("Chromimum_{cellranger_organism}_Transcriptome_Probe_Set_v1.1.0_GRCh38-2024-A.csv")
      ),
      "10xflex_v1.1_multi" = file.path(
        flex_probe_dir,
        glue::glue("Chromimum_{cellranger_organism}_Transcriptome_Probe_Set_v1.1.0_GRCh38-2024-A.csv")
      )
    )
  }

  # add directory for visium probes
  if (include_visium) {
    # override the assembly for naming these files
    visium_assembly <- assembly
    if (assembly == "GRCm38") {
      visium_assembly <- "mm10"
    }

    # determine which visium probes to include
    include_techs <- visium_probe_map |>
      dplyr::filter(
        organism == cellranger_organism,
        internal_assembly == assembly,
        ensembl_version == version
      ) |>
      dplyr::pull(technology)

    # create all but keep only the include_techs
    visium_probe_files <- list(
      "visium1_v1" = file.path(
        visium_probe_dir,
        glue::glue("Visium_{cellranger_organism}_Transcriptome_Probe_Set_v1.0_{visium_assembly}-2020-A.csv")
      ),
      "visium2_v2" = file.path(
        visium_probe_dir,
        glue::glue("Visium_{cellranger_organism}_Transcriptome_Probe_Set_v2.0_{visium_assembly}-2020-A.csv")
      ),
      "visium2_v2.1" = file.path(
        visium_probe_dir,
        glue::glue("Visium_{cellranger_organism}_Transcriptome_Probe_Set_v2.1.0_{visium_assembly}-2024-A.csv")
      ),
      "visium-hd_v2" = file.path(
        visium_probe_dir,
        glue::glue("Visium_{cellranger_organism}_Transcriptome_Probe_Set_v2.0_{visium_assembly}-2020-A.csv")
      ),
      "visium-hd_v2.1" = file.path(
        visium_probe_dir,
        glue::glue("Visium_{cellranger_organism}_Transcriptome_Probe_Set_v2.1_{visium_assembly}-2024-A.csv")
      )
    ) |>
      purrr::keep_at(include_techs)

    # add to json
    if (length(visium_probe_files) > 0) {
      json_entry$visium_probe_files <- visium_probe_files
    } else {
      json_entry$visium_probe_files <- ""
    }
  }
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
