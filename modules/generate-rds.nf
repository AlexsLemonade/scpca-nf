
// generate unfiltered and filtered RDS file using scpcaTools
process generate_rds{
    container 'ghcr.io/alexslemonade/scpca-r'
    publishDir "${params.outdir}/${meta.sample_id}/${meta.library_id}"
    input: 
        tuple val(meta), path(alevin_dir)
    output:
        tuple val(meta), path(unfiltered_rds), path(filtered_rds)
    script:
    unfiltered_rds = "${params.outdir}/${meta.sample_id}/${meta.library_id}/unfiltered.rds"
    filtered_rds = "${params.outdir}/${meta.sample_id}/${meta.library_id}/filtered.rds"
    """
    Rscript --vanilla scripts/generate_output_files.R \
      -s ${meta.seq_unit} \
      -a ${alevin_dir} \
      -u ${unfiltered_rds} \
      -f ${filtered_rds}
    """
}