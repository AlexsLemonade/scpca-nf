#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// integration specific parameters
params.integration_metafile = 's3://ccdl-scpca-data/sample_info/scpca-integration-metadata.tsv'
params.integration_group = "All"

// define path to integration template
integration_template = "${projectDir}/templates/integration-report.Rmd"

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
  label 'mem_16'
  publishDir "${params.checkpoints_dir}/merged_sces"
  input:
    tuple val(integration_group), val(library_ids), path(scpca_nf_file)
  output:
    tuple val(integration_group), path(merged_sce_file)
  script:
    input_library_ids = library_ids.join(',')
    input_sces = scpca_nf_file.join(',')
    merged_sce_file = "${integration_group}_merged.rds"
    """
    merge_sces.R \
      --input_library_ids "${input_library_ids}" \
      --input_sce_files "${input_sces}" \
      --output_sce_file "${merged_sce_file}" \
      --n_hvg ${params.num_hvg} \
      --threads ${task.cpus}
    """
  stub:
    merged_sce_file = "${integration_group}_merged.rds"
    """
    touch ${merged_sce_file}
    """

}

// integrate with fastMNN
process integrate_fastmnn {
  container params.SCPCATOOLS_CONTAINER
  label 'mem_16'
  label 'cpus_4'
  input:
    tuple val(integration_group), path(merged_sce_file)
  output:
    tuple val(integration_group), path(integrated_sce_file)
  script:
    integrated_sce_file = "${integration_group}.rds"
    """
    integrate_sce.R \
      --input_sce_file "${merged_sce_file}" \
      --output_sce_file "${integrated_sce_file}" \
      --method "fastMNN" \
      --seed ${params.seed} \
      --threads ${task.cpus}
    """
  stub:
    integrated_sce_file = "${integration_group}.rds"
    """
    touch ${integrated_sce_file}
    """
}

// integrate with fastMNN
process integrate_harmony {
  container params.SCPCATOOLS_CONTAINER
  publishDir "${params.results_dir}/integration/${integration_group}"
  label 'mem_16'
  input:
    tuple val(integration_group), path(merged_sce_file)
  output:
    tuple val(integration_group), path(integrated_sce_file)
  script:
    integrated_sce_file = "${integration_group}.rds"
    """
    integrate_sce.R \
      --input_sce_file "${merged_sce_file}" \
      --output_sce_file "${integrated_sce_file}" \
      --method "harmony" \
      --seed ${params.seed}
    """
  stub:
    integrated_sce_file = "${integration_group}.rds"
    """
    touch ${integrated_sce_file}
    """
}

// create integrated report and single object
process integration_report {
  container params.SCPCATOOLS_CONTAINER
  publishDir "${params.results_dir}/integration/${integration_group}"
  label 'mem_16'
  input:
    tuple val(integration_group), path(integrated_sce_file)
    path(report_template)
  output:
    path(integration_report)
  script:
    integration_report = "${integration_group}_summary_report.html"
    """
    Rscript -e "rmarkdown::render('${report_template}', \
                                  output_file = '${integration_report}', \
                                  params = list(integration_group = '${integration_group}', \
                                                integrated_sce = '${integrated_sce_file}', \
                                                batch_column = 'library_id'))"
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
    libraries_ch = Channel.fromPath(params.run_metafile)
      .splitCsv(header: true, sep: '\t')
      // only include single-cell/single-nuclei and make sure no CITE-seq/ hashing libraries
      .filter{it.seq_unit in ['cell', 'nucleus']}
      .map{[
        library_id: it.scpca_library_id,
        scpca_nf_file: "${params.results_dir}/${it.scpca_project_id}/${it.scpca_sample_id}/${it.scpca_library_id}_processed.rds"
      ]}
      .unique()

    grouped_meta_ch = integration_meta_ch
      .map{[it.library_id, it.integration_group]}
      // pull out library_id from meta and use to join
      .combine(libraries_ch.map{[it.library_id, it.scpca_nf_file]}, by: 0)
      // create tuple of integration group, library ID, and output file from scpca_nf
      .map{[
        it[1], // integration_group
        it[0], // library_id
        file(it[2]) // scpca_nf_file
        ]}
      // grouped tuple of [integration_group, [library_id1, library_id2, ...], [sce_file1, sce_file2, ...]]
      .groupTuple(by: 0)

    merge_sce(grouped_meta_ch)

    // integrate using fastmnn
    integrate_fastmnn(merge_sce.out)

    // integrate using harmony
    integrate_harmony(integrate_fastmnn.out)

    // generate integration report
    integration_report(integrate_harmony.out, file(integration_template))
}
