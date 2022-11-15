
// generate QC report from unfiltered and filtered SCE.rds files using scpcaTools

process sce_qc_report{
    container params.SCPCATOOLS_CONTAINER
    tag "${meta.library_id}"
    publishDir "${params.results_dir}/${meta.project_id}/${meta.sample_id}", mode: 'copy'
    input:
        tuple val(meta), path(unfiltered_rds), path(filtered_rds), path(processed_rds)
    output:
        tuple val(meta), path(qc_report), path(metadata_json)
    script:
        qc_report = "${meta.library_id}_qc.html"
        metadata_json = "${meta.library_id}_metadata.json"
        workflow_url = workflow.repository ?: workflow.manifest.homePage
        workflow_version = workflow.revision ?: workflow.manifest.version
        """
        sce_qc_report.R \
          --library_id "${meta.library_id}" \
          --sample_id "${meta.sample_id}" \
          --unfiltered_sce ${unfiltered_rds} \
          --filtered_sce ${filtered_rds} \
          --qc_report_file ${qc_report} \
          --metadata_json ${metadata_json} \
          --technology "${meta.technology}" \
          --seq_unit "${meta.seq_unit}" \
          --genome_assembly "${meta.ref_assembly}" \
          --workflow_url "${workflow_url}" \
          --workflow_version "${workflow_version}" \
          --workflow_commit "${workflow.commitId}"
        """
}
