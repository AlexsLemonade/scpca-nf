#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// integration specific parameters
params.integration_metafile = 's3://ccdl-scpca-data/sample_info/scpca-integration-metadata.tsv'
params.integration_group = "All"

// parameter checks
param_error = false

if (!file(params.run_metafile).exists()) {
  log.error("The 'run_metafile' file '${params.run_metafile}' can not be found.")
  param_error = true
}

if (!file(params.integration_metafile).exists()) {
  log.error("The 'integration_metafile' file '${params.integration_metafile}' can not be found.")
  param_error = true
}

if(param_error){
  System.exit(1)
}

// merge individual SCE objects into one SCE object
process merge_sce {
  container params.SCPCATOOLS_CONTAINER
  publishDir "${params.checkpoints_dir}/merged_sces"
  input:
    tuple val(integration_group), val(library_ids), path(scpca_nf_file)
  output:
    path merged_sce_file
  script:
    input_sces = scpca_nf_file.join(',')
    input_library_ids = library_ids.join(',')
    merged_sce_file = "${integration_group}_merged.txt"
    """
    echo $library_ids $input_sces > $merged_sce_file
    # would then add in a call to merging script with the text file as input
    # the merging script would output the merged_sce_file

    """

}

workflow {

    // select projects to integrate from params
    integration_groups = params.integration_group?.tokenize(',') ?: []
    integration_groups_all = integration_groups[0] == "All" // create logical for including all groups or not when filtering later

    // create channel of integration group and libraries to integrate
    integration_meta_ch = Channel.fromPath(params.integration_metafile)
      .splitCsv(header: true, sep: '\t')
      .map{[
        library_id: it.scpca_library_id,
        integration_group: it.integration_group,
        submitter: it.submitter
      ]}
      .filter{integration_groups_all  || (it.integration_group in integration_groups)}

    // channel with run metadata, keeping only the columns we need
    runs_ch = Channel.fromPath(params.run_metafile)
      .splitCsv(header: true, sep: '\t')
      .map{[
        library_id: it.scpca_library_id,
        scpca_nf_file: "${params.results_dir}/${it.scpca_sample_id}/${it.scpca_library_id}_processed.rds"
      ]}

    all_meta_ch = integration_meta_ch
      .map{[it["library_id"]] + it }
      // pull out library_id from meta and use to join
      .combine(runs_ch.map{[it["library_id"]] + it }, by: 0)
      // create tuple of integration group, library ID, and output file from scpca_nf
      .map{[
        it[1].integration_group,
        it[0], // library_id
        file(it[2].scpca_nf_file)
        ]}
      // grouped tuple of [integration_group, [file1, file2, file3, ...]]
      .groupTuple(by: 0)

    merge_sce(all_meta_ch)

}

