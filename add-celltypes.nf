#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { annotate_celltypes } from './modules/classify-celltypes.nf'

// parameter checks
param_error = false

if (!file(params.run_metafile).exists()) {
  log.error("The 'run_metafile' file '${params.run_metafile}' can not be found.")
  param_error = true
}

if (!file(params.celltype_refs_metafile).exists()) {
  log.error("The 'celltype_refs_metafile' file '${params.celltype_refs_metafile}' can not be found.")
  param_error = true
}

workflow {

  // select runs to use
  if (params.project){
    // projects will use all runs in the project & supersede run_ids
    run_ids = []
  }else{
    run_ids = params.run_ids?.tokenize(',') ?: []
  }
  run_all = run_ids[0] == "All"
  if (run_all){
    log.info("Executing workflow for all runs in the run metafile.")
  }

  // read in metadata file and filter to libraries/ projects of interest
  processed_sce_ch = Channel.fromPath(params.run_metafile)
    .splitCsv(header: true, sep: '\t')
    .map{[
        run_id: it.scpca_run_id,
        library_id: it.scpca_library_id,
        sample_id: it.scpca_sample_id,
        project_id: it.scpca_project_id,
        submitter: it.submitter,
        technology: it.technology,
        seq_unit: it.seq_unit,
    ]}
    .filter{it.seq_unit in ['cell', 'nucleus']}
    // filter to only single-cell and remove any CITE-seq or multiplexed data
    .filter{it.technology.startsWith("10Xv")}
    .filter{run_all
             || (it.run_id in run_ids)
             || (it.library_id in run_ids)
             || (it.sample_id in run_ids)
             || (it.submitter == params.project)
             || (it.project_id == params.project)
            }
    // tuple of meta, processed rds file to use as input to cell type annotation
    .map{meta -> tuple(meta,
                       file("${params.results_dir}/${meta.project_id}/${meta.sample_id}/${meta.library_id}_processed.rds")
                       )}

    annotate_celltypes(processed_sce_ch)
}
