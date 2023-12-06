#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Workflow to merge SCE objects into a single object.
// This workflow does NOT perform integration, i.e. batch correction.

// define path to merge template
// TODO: Establish this merge-report.Rmd file
//merge_template = "${projectDir}/templates/merge-report.Rmd"

// parameter checks
param_error = false

// check that at least one project has been provided
if(!params.project) {
  log.error("At least one 'project' must be specified for merging.")
  param_error = true
}

// check for provided run file
if (!file(params.run_metafile).exists()) {
  log.error("The 'run_metafile' file '${params.run_metafile}' can not be found.")
  param_error = true
}

if(param_error){
  System.exit(1)
}

// merge individual SCE objects into one SCE object
process merge_sce {
  container params.SCPCATOOLS_CONTAINER
  label 'mem_16'
  publishDir "${params.checkpoints_dir}/merged"
  input:
    tuple val(project_id), val(library_ids), path(scpca_nf_file)
  output:
    tuple val(project_id), path(merged_sce_file)
  script:
    input_library_ids = library_ids.join(',')
    input_sces = scpca_nf_file.join(',')
    merged_sce_file = "${project_id}_merged.rds"
    """
    merge_sces.R \
      --input_library_ids "${input_library_ids}" \
      --input_sce_files "${input_sces}" \
      --output_sce_file "${merged_sce_file}" \
      --n_hvg ${params.num_hvg} \
      --threads ${task.cpus}
    """
  stub:
    merged_sce_file = "${merge_group}_merged.rds"
    """
    touch ${merged_sce_file}
    """

}

// create merge report
process merge_report {
  container params.SCPCATOOLS_CONTAINER
  publishDir "${params.results_dir}/merged/${merge_group}"
  label 'mem_16'
  input:
    tuple val(merge_group), path(merged_sce_file)
    path(report_template)
  output:
    path(merge_report)
  script:
    merge_report = "${merge_group}_summary_report.html"
    """
    Rscript -e "rmarkdown::render( \
      '${report_template}', \
      output_file = '${merge_report}', \
      params = list(merge_group = '${merge_group}', \
                    merged_sce = '${merged_sce_file}', \
                    batch_column = 'library_id') \
      )"
    """
  stub:
    merge_report = "${merge_group}_summary_report.html"
    """
    touch ${merge_report}
    """
}

workflow {

    // grab project ids to run
    project_ids = params.project?.tokenize(',') ?: []

    // read in run metafile, filter to projects of interest, and group by project
    grouped_libraries_ch = Channel.fromPath(params.run_metafile)
      .splitCsv(header: true, sep: '\t')
      // filter to only include specified project ids
      .filter{it.scpca_project_id in project_ids}
      // only include single-cell/single-nuclei and make sure no CITE-seq/ hashing libraries
      .filter{it.seq_unit in ['cell', 'nucleus']}
      // create tuple of [project id, library_id, processed_sce_file]
      .map{[
        it.scpca_project_id,
        it.scpca_library_id,
        "${params.results_dir}/${it.scpca_project_id}/${it.scpca_sample_id}/${it.scpca_library_id}_processed.rds"
      ]}
      // only include samples that have been processed through scpca-nf
      .filter{file(it[2]).exists()}
      // make sure we don't have any duplicates of the same library ID hanging around
      // this shouldn't be the case since we removed CITE-seq and cell-hashing
      .unique()
      // group tuple by project id, [project_id, [library_id1, library_id2, ...], [sce_file1, sce_file2, ...]]
      .groupTuple(by: 0)

    merge_sce(grouped_libraries_ch)

    // TODO: generate merge report
    //merge_report(merge_sce.out, file(merge_template))
}
