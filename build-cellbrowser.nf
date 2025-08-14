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
  if (!params.cellbrowser_dir) {
    log.error("The 'cellbrowser_dir' directory is required for generating cellbrowser output")
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
  //
  libraries_ch = Channel.fromPath(params.run_metafile)
    .splitCsv(header: true, sep: '\t')
    // filter to run all ids or just specified ones
    .map{it -> [
        project_id: it.scpca_project_id,
        library_id: it.scpca_library_id,
        sample_id: it.scpca_sample_id.split(";").sort().join(","),
        run_id: it.scpca_run_id
    ]}
    .filter{
      run_all
      || (it.run_id in run_ids)
      || (it.library_id in run_ids)
      || (it.sample_id in run_ids)
      || (it.project_id in project_ids)
    }
    .unique{ it.library_id }
    .map{it -> [
      it, // meta
      file("${params.results_dir}/${it.project_id}/${it.sample_id}/${it.library_id}_processed_rna.h5ad")
    ]}
    // only include libraries where the h5ad file exists
    .filter{ it[1].exists() }

  cellbrowser_build(libraries_ch)
}
