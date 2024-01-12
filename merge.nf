#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Workflow to merge SCE objects into a single object.
// This workflow does NOT perform integration, i.e. batch correction.

// define path to merge template
merge_template = "${projectDir}/templates/merge-report.rmd"

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
  label 'mem_32'
  publishDir "${params.results_dir}/merged/${merge_group_id}"
  input:
    tuple val(merge_group_id), val(has_adt), val(multiplexed), val(library_ids), path(scpca_nf_file)
  output:
    tuple path(merged_sce_file), val(merge_group_id), val(has_adt), val(multiplexed)
  script:
    input_library_ids = library_ids.join(',')
    input_sces = scpca_nf_file.join(',')
    merged_sce_file = "${merge_group_id}_merged.rds"
    """
    merge_sces.R \
      --input_library_ids "${input_library_ids}" \
      --input_sce_files "${input_sces}" \
      --output_sce_file "${merged_sce_file}" \
      --n_hvg ${params.num_hvg} \
      ${has_adt ? "--include_altexp" : ''} \
      ${multiplexed ? "--multiplexed" : '' } \
      --threads ${task.cpus}
    """
  stub:
    merged_sce_file = "${merge_group_id}_merged.rds"
    """
    touch ${merged_sce_file}
    """

}

// create merge report
process generate_merge_report {
  container params.SCPCATOOLS_CONTAINER
  publishDir "${params.results_dir}/merged/${merge_group_id}"
  label 'mem_16'
  input:
    tuple path(merged_sce_file), val(merge_group_id), val(has_adt), val(multiplexed)
    path(report_template)
  output:
    path(merge_report)
  script:
    merge_report = "${merge_group_id}_summary_report.html"
    """
    Rscript -e "rmarkdown::render( \
      '${report_template}', \
      output_file = '${merge_report}', \
      params = list(merge_group = '${merge_group_id}', \
                    merged_sce_file = '${merged_sce_file}', \
                    batch_column = 'library_id') \
      )"
    """
  stub:
    merge_report = "${merge_group_id}_summary_report.html"
    """
    touch ${merge_report}
    """
}

process export_anndata{
    container params.SCPCATOOLS_CONTAINER
    label 'mem_32'
    tag "${merge_group_id}"
    publishDir "${params.results_dir}/merged/${merge_group_id}", mode: 'copy'
    input:
      tuple path(merged_sce_file), val(merge_group_id), val(has_adt)
    output:
      tuple val(merge_group_id), path("${merge_group_id}_merged_*.hdf5")
    script:
      rna_hdf5_file = "${merge_group_id}_merged_rna.hdf5"
      feature_hdf5_file = "${merge_group_id}_merged_adt.hdf5"
      """
      sce_to_anndata.R \
        --input_sce_file ${merged_sce_file} \
        --output_rna_h5 ${rna_hdf5_file} \
        --output_feature_h5 ${feature_hdf5_file} \
        ${has_adt ? "--feature_name adt" : ''}

      # move normalized counts to X in AnnData
      move_counts_anndata.py --anndata_file ${rna_hdf5_file}
      ${has_adt ? "move_counts_anndata.py --anndata_file ${feature_hdf5_file}" : ''}
      """
    stub:
      rna_hdf5_file = "${merge_group_id}_merged_rna.hdf5"
      feature_hdf5_file = "${merge_group_id}_merged_adt.hdf5"
      """
      touch ${rna_hdf5_file}
      ${has_adt ? "touch ${feature_hdf5_file}" : ''}
      """
}

workflow {

    // grab project ids to run
    project_ids = params.project?.tokenize(',') ?: []

    // read in run metafile and filter to projects of interest
    libraries_ch = Channel.fromPath(params.run_metafile)
      .splitCsv(header: true, sep: '\t')
      // filter to only include specified project ids
      .filter{it.scpca_project_id in project_ids}

    // get all projects that contain at least one library with CITEseq
    adt_projects = libraries_ch
      .filter{it.technology.startsWith('CITEseq')}
      .collect{it.scpca_project_id}
      .unique()

    multiplex_projects = libraries_ch
      .filter{it.technology.startsWith('cellhash')}
      .collect{it.scpca_project_id}
      .unique()

    grouped_libraries_ch = libraries_ch
      // only include single-cell/single-nuclei which ensures we don't try to merge libraries from spatial or bulk data
      .filter{it.seq_unit in ['cell', 'nucleus']}
      // create tuple of [project id, library_id, processed_sce_file]
      .map{[
        it.scpca_project_id,
        it.scpca_library_id,
        file("${params.results_dir}/${it.scpca_project_id}/${it.scpca_sample_id}/${it.scpca_library_id}_processed.rds")
      ]}
      // only include libraries that have been processed through scpca-nf
      .filter{file(it[2]).exists()}
      // only one row per library ID, this removes all the duplicates that may be present due to CITE/hashing
      .unique()
      // group tuple by project id: [project_id, [library_id1, library_id2, ...], [sce_file1, sce_file2, ...]]
      .groupTuple(by: 0)
      // add in boolean for if project contains samples with adt
      .map{project_id, library_id_list, sce_file_list -> tuple(
        project_id,
        project_id in adt_projects, // determines if altExp should be included in the merged object
        project_id in multiplex_projects, // determines if sample metadata should be added to colData and to skip anndata
        library_id_list,
        sce_file_list
      )}

    merge_sce(grouped_libraries_ch)

    // generate merge report
    generate_merge_report(merge_sce.out, file(merge_template))

    // export merged objects to AnnData
    anndata_ch = merge_sce.out
      .filter{!it[3]} // remove multiplexed samples before export
      .map{it.take(3)} // keep everything but multiplexed logical


    export_anndata(anndata_ch)
}
