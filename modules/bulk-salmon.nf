#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process fastp{
    container params.FASTP_CONTAINER
    label 'cpus_8'
    tag "${meta.run_id}-bulk"
    publishDir "${params.outdir}/internal/fastp/${meta.library_id}"
    input: 
        tuple val(meta),path(read1), path(read2)
    output: 
        tuple val(meta),path(trimmed_read1),path(trimmed_read2)
    script: 
        trimmed_read1 = "${meta.run_id}-trimmed-R1.fastq.gz"
        trimmed_read2 = "${meta.run_id}-trimmed-R2.fastq.gz"
        """
        fastp --in1 ${read1} \
        ${meta.library == 'paired_end' ? "--in2 ${read2}":""} \
        --out1 "${trimmed_read1} \
        ${meta.library == 'paired_end' ? "--out2 ${trimmed_read2}":""} \
        --html "${id}_fastp.html" \
        --json "${id}_fastp.json" \
        --trim_poly_g \
        --report_title ${meta.library_id}
        """

}

process salmon{
    container params.SALMON_CONTAINER
    label 'cpus_8'
    tag "${meta.run_id}-bulk"
    publishDir "${params.outdir}/internal/salmon/${meta.library_id}"
    input: 
        tuple val(meta),path(trimmed_read1),path(trimmed_read2)
    output: 
        tuple val(meta),path(salmon_results)
    script:
        salmon_results = "${meta.run_id}-salmon"
        """
        salmon quant -i ${index} /
        -l A /
        -1 ${trimmed_read1} /
        ${meta.library == 'paired_end' ? "-2 ${trimmed_read2}":""} /
        -o ${salmon_results} /
        --threads ${task.cpus}
        """

}

workflow bulk_quant_rna {
    take: bulk_channel 
    // a channel with a map of metadata for each rna library to process
    main: 
    // create tuple of (metadata map, [Read 1 files], [Read 2 files])
    bulk_reads_ch = bulk_channel
      .map{meta -> tuple(meta,
                         file("s3://${meta.s3_prefix}/*_R1_*.fastq.gz"),
                         file("s3://${meta.s3_prefix}/*_R2_*.fastq.gz"))}

    fastp(bulk_reads_ch) \
      | salmon
    
    emit: salmon.out
}
