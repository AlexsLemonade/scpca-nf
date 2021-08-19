
// generate unfiltered and filtered RDS file using scpcaTools
process generate_rds{
    container 'ghcr.io/alexslemonade/scpca-r'
    label 'cpus_8'
    publishDir "${params.outdir}/${meta.sample_id}"
    input: 
        tuple val(meta), path(alevin_dir)
    output:
        tuple val(meta), path(unfiltered_rds), path(filtered_rds)
    script:
        unfiltered_rds = "${meta.library_id}_unfiltered.rds"
        filtered_rds = "${meta.library_id}_filtered.rds"
        """
        generate_output_files.sh \
          --seq_unit ${meta.seq_unit} \
          --alevin_dir ${alevin_dir} \
          --unfiltered_file ${unfiltered_rds} \
          --filtered_file ${filtered_rds}
        """
}
