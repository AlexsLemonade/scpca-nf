
// generate QC report from unfiltered and filtered SCE.rds files using scpcaTools

process sce_qc_report{
    container params.SCPCATOOLS_CONTAINER
    publishDir "${params.outdir}/publish/${meta.project_id}/${meta.sample_id}"
    input: 
        tuple val(meta), path(unfiltered_rds), path(filtered_rds)
    output:
        tuple val(meta), path(qc_report)
    script:
        qc_report = "${meta.library_id}_qc.html"
        """
        sce_qc_report.R \
          --sample_id {meta.library_id} \
          --unfiltered_sce ${unfiltered_rds} \
          --filtered_sce ${filtered_rds} \
          --output_file ${qc_report}
        """
}
