#!/usr/bin/env nextflow
nextflow.enable.dsl=2

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

process merge_sce {
  container params.SCPCATOOLS_CONTAINER
  publishDir "${params.checkpoints_dir}/merged_sces"
  input:
    tuple val(int_group), path(scpca_nf_file), val(library_ids), val(sample_ids)
  output:
    path merged_sce_file
  script:
    merged_sce_file = "${int_group}_merged.rds"
    """
    echo $scpca_nf_file > $merged_sce_file
    """

}

workflow {

    // select projects to integrate from params
    int_groups = params.integration_group?.tokenize(',') ?: []
    int_groups_all = int_groups[0] == "All" // create logical for including all groups or not when filtering later

    // create channel of integration group and libraries to integrate
    integration_meta_ch = Channel.fromPath(params.integration_metafile)
      .splitCsv(header: true, sep: '\t')
      .map{[
        library_id: it.scpca_library_id,
        integration_group: it.integration_group,
        submitter: it.submitter
      ]}
      .filter{int_groups_all  || (it.integration_group in int_groups)}

    // channel with run metadata, keeping only the columns we need
    runs_ch = Channel.fromPath(params.run_metafile)
      .splitCsv(header: true, sep: '\t')
      .map{[
        library_id: it.scpca_library_id,
        run_id: it.scpca_run_id,
        sample_id: it.scpca_sample_id.split(";").sort().join(","),
        scpca_nf_file: "${params.outdir}/${it.scpca_sample_id}/${it.scpca_library_id}_processed.rds"
      ]}
      .map{meta -> tuple(meta,
                         scpca_nf_file: file(meta.scpca_nf_file)
                         )}

    all_meta_ch = integration_meta_ch
      .map{[it["library_id"]] + it }
      // pull out library_id from meta and use to join
      .join(runs_ch.map{[it[1]["library_id"]] + it }, by: 0)
      // create tuple of just integration group and output file from scpca_nf
      .map{[
        it[1].integration_group,
        it[2].scpca_nf_file,
        it[3].library_id,
        it[3].sample_id]}
      // grouped tuple of [integration_group, [file1, file2, file3, ...],
      //                    [library1, library2, library3...], [sample1, sample2, sample3...]]
      .groupTuple(by: 0)

    merge_sce(all_meta_ch)

}

