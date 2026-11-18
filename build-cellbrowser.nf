#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { cellbrowser_build } from './modules/cellbrowser.nf'

def check_parameters() {
  // parameter check function
  def param_error = false
  if (!file(params.run_metafile).exists()) {
    log.error("The 'run_metafile' file '${params.run_metafile}' can not be found.")
    param_error = true
  }

  if (!file(params.sample_metafile).exists()) {
    log.error("The 'sample_metafile' file '${params.sample_metafile}' can not be found.")
    param_error = true
  }

  // cell browser config checks
  if (!params.cellbrowser_dirname) {
    log.error("The 'cellbrowser_dirname' directory is required for generating cellbrowser output")
    param_error = true
  }
  if (!file(params.project_metafile).exists()) {
    log.error("The 'project_metafile' file '${params.project_metafile}' can not be found.")
    param_error = true
  }

  if (param_error) {
    System.exit(1)
  }
}

workflow {
  check_parameters()

  // select runs to use
  def run_ids = []
  def project_ids = []
  if (params.project) {
    // projects will use all runs in the project & supersede run_ids
    // allow for processing of multiple projects at once
    project_ids = params.project?.tokenize(',') ?: []
  } else {
    run_ids = params.run_ids?.tokenize(',') ?: []
  }
  def run_all = run_ids[0] == "All"
  if (!run_all) {
    log.warn("Some runs are not being included in the Cell Browser output; is this correct?")
  }
  // technologies with multiple samples per library, which are split into one entry per library-sample combination
  def flex_technologies = ["10xflex_v1.1_multi"]

  libraries_ch = channel.fromPath(params.run_metafile)
    .splitCsv(header: true, sep: '\t')
    // split each row into one or more library-sample combinations depending on the technology
    .flatMap{ row ->
        def technology = row.technology.toLowerCase()
        def sample_ids = row.scpca_sample_id.split("[;,]").sort()
        def library_meta = [
          project_id: row.scpca_project_id,
          run_id: row.scpca_run_id,
          library_id: row.scpca_library_id,
          technology: technology
        ]

        if (technology in flex_technologies) {
          // flex: one row per library-sample pair
          return sample_ids.collect{ sample_id ->
            library_meta + [
              sample_id: sample_id,
              unique_id: "${row.scpca_library_id}-${sample_id}".toString()
            ]
          }
        } else {
          // everything else, including cellhash multiplexed libraries: one row per library
          return [
            library_meta + [
              sample_id: sample_ids.join(","),
              unique_id: row.scpca_library_id
            ]
          ]
        }
      }
    // filter after splitting flex libraries so individual sample ids can be matched
    .filter{ it ->
      run_all
      || (it.run_id in run_ids)
      || (it.library_id in run_ids)
      || it.sample_id.split(",").any{ sample_id -> sample_id in run_ids }
      || (it.project_id in project_ids)
    }
    .unique{ it.unique_id }
    .map{ meta ->
      def h5ad_file = file("${params.results_dir}/${meta.project_id}/${meta.sample_id}/${meta.unique_id}_processed_rna.h5ad")
      [meta, h5ad_file]
    }
    // only include libraries where the h5ad file exists
    .filter{ it[1].exists() }

  cellbrowser_build(libraries_ch)
}
