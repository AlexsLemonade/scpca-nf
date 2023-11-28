#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Workflow to merge SCE objects into a single object.
// This workflow does NOT perform integration, i.e. batch correction.


// merge-specific parameters
// TODO: Update approach to define merge groupings.
params.merge_metafile = 's3://ccdl-scpca-data/sample_info/scpca-integration-metadata.tsv'
params.merge_group = "All"

// define path to merge template
// TODO: Establish this merge-report.Rmd file
//merge_template = "${projectDir}/templates/merge-report.Rmd"

// parameter checks
param_error = false

if (!file(params.run_metafile).exists()) {
  log.error("The 'run_metafile' file '${params.run_metafile}' can not be found.")
  param_error = true
}

if (!file(params.merge_metafile).exists()) {
  log.error("The 'merge_metafile' file '${params.merge_metafile}' can not be found.")
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
    tuple val(merge_group), val(library_ids), path(scpca_nf_file)
  output:
    tuple val(merge_group), path(merged_sce_file)
  script:
    input_library_ids = library_ids.join(',')
    input_sces = scpca_nf_file.join(',')
    merged_sce_file = "${merge_group}_merged.rds"
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

// create merge report and single object
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

}

workflow {

    // select projects to merge from params
    merge_groups = params.merge_group?.tokenize(',') ?: []
    merge_groups_all = merge_groups[0] == "All" // create logical for including all groups or not when filtering later

    // create channel of merge group and libraries to merge
    merge_meta_ch = Channel.fromPath(params.merge_metafile)
      .splitCsv(header: true, sep: '\t')
      .map{[
        library_id: it.scpca_library_id,
        merge_group: it.merge_group,
        submitter: it.submitter
      ]}
      .filter{merge_groups_all  || (it.merge_group in merge_groups)}

    // channel with run metadata, keeping only the columns we need
    libraries_ch = Channel.fromPath(params.run_metafile)
      .splitCsv(header: true, sep: '\t')
      // only include single-cell/single-nuclei and make sure no CITE-seq/ hashing libraries
      .filter{it.seq_unit in ['cell', 'nucleus']}
      .map{[
        library_id: it.scpca_library_id,
        scpca_nf_file: "${params.results_dir}/${it.scpca_project_id}/${it.scpca_sample_id}/${it.scpca_library_id}_processed.rds"
      ]}
      .unique()

    grouped_meta_ch = merge_meta_ch
      .map{[it.library_id, it.merge_group]}
      // pull out library_id from meta and use to join
      .combine(libraries_ch.map{[it.library_id, it.scpca_nf_file]}, by: 0)
      // create tuple of merge group, library ID, and output file from scpca_nf
      .map{[
        it[1], // merge_group
        it[0], // library_id
        file(it[2]) // scpca_nf_file
        ]}
      // grouped tuple of [merge_group, [library_id1, library_id2, ...], [sce_file1, sce_file2, ...]]
      .groupTuple(by: 0)

    merge_sce(grouped_meta_ch)

    // TODO: generate merge report
    //merge_report(merge_sce.out, file(merge_template))
}
