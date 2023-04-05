
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

if(!file.exists(opt$ref_metadata)){
  stop("ref_metadata file does not exist.")
}

ref_metadata <- readr::read_tsv(opt$ref_metadata) |>
  dplyr::mutate(reference_name = glue::glue("{organism}.{assembly}.{version}"))

# check for columns

ref_rootdir <- file.path("s3://scpca-references")

create_ref_entry <- function(organism,
                             assembly,
                             version,
                             reference_name,
                             ref_rootdir){

  ref_dir <- file.path(ref_rootdir,
                       tolower(organism),
                       glue::glue("ensembl-{version}"))

  json_entry <- list(
    reference_name = reference_name,
    ref_fasta = file.path(ref_dir, "fasta",
                           glue::glue("{organism}.{assembly}.dna.primary_assembly.fa.gz")),
    ref_fasta_index = file.path(ref_dir, "fasta",
                                 glue::glue("{organism}.{assembly}.dna.primary_assembly.fa.fai")),
    ref_gtf = file.path(ref_dir, "annotation",
                        glue::glue("{reference_name}.gtf.gz")),
    splici_index = file.path(ref_dir, "salmon_index",
                             glue::glue("{reference_name}.spliced_intron.txome")),
    bulk_index = file.path(ref_dir, "salmon_index",
                            glue::glue("{reference_name}.spliced_cdna.txome")),
    cellranger_index = file.path(ref_dir, "cellranger_index",
                                  glue::glue("{reference_name}_cellranger_full")),
    star_index = file.path(ref_dir, "star_index",
                            glue::glue("{reference_name}.star_idx")),
    mito_file = file.path(ref_dir, "annotation",
                           glue::glue("{reference_name}.mitogenes.txt")),
    t2g_3col_path = file.path(ref_dir, "annotation",
                               glue::glue("{reference_name}.spliced_intron.tx2gene_3col.tsv")),
    t2g_bulk_path = file.path(ref_dir, "annotation",
                               glue::glue("{reference_name}.spliced_cdna.tx2gene.tsv"))
    ) |>
    purrr::map(jsonlite::unbox)

    return(json_entry)
}

all_entries <- purrr::pmap(list(ref_metadata$organism,
                                ref_metadata$assembly,
                                ref_metadata$version,
                                ref_metadata$reference_name),
                           \(organism, assembly, version, reference_name)
                           create_ref_entry(organism, assembly, version, reference_name, ref_rootdir)) |>
  purrr::set_names(ref_metadata$reference_name)

# Write to JSON
jsonlite::write_json(all_entries,
                     path = opt$ref_json,
                     pretty = TRUE)
