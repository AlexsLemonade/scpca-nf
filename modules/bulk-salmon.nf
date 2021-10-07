#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process fastp{
    container params.FASTP_CONTAINER
    label 'cpus_8'
    tag "${meta.run_id}-bulk"
    publishDir "${params.outdir}/internal/fastp/${meta.library_id}"
    input: 
        tuple val(meta),path(read1),path(read2)
    output: 
        tuple val(meta),path(trimmed_read1),path(trimmed_read2)
    script: 
        trimmed_data = "${params.outdir}/data/trimmed"
        fastp_reports = "${params.outdir}/reports/fastp"
        """
        mkdir -p data/trimmed
        mkdir -p reports/fastp

        fastp --in1 {input.r1} --in2 {input.r2} /
        --out1 "${trimmed_data}/${id}_R1_001_fastq.gz" /
        --out2 "${trimmed_data}/${id}_R2_001_fastq.gz" /
        --html "${fastp_reports}/${id}_fastp.html" /
        --json "${fastp_reports}/${id}_fastp.json" /
        --trim_poly_g /
        --report_title '${id} report'
        """

}

process salmon{
    container 'quay.io/biocontainers/salmon:1.4.0--hf69c8f4_0'
    label 'cpus_8'
    tag "${id}-${index}"
    publishDir "${params.outdir}"
    input: 
        tuple val(meta),path(trimmed_read1),path(trimmed_read2)
    output: 
        tuple val(meta),path()
    script:
        salmon_results = "${params.outdir}/salmon-quant"
        """
        salmon quant -i ${index} /
        -l A /
        -1 "${trimmed_data}/${id}_R1_001_fastq.gz" /
        -2 "${trimmed_data}/${id}_R2_001_fastq.gz" /
        -o ${salmon_results} /
        --threads ${task.cpus}
        """

}

workflow bulk_quant_rna {
    take: bulk_channel 
    // a channel with a map of metadata for each rna library to process
    main: 
    // create tuple of (metadata map, [Read 1 files], [Read 2 files])
    single_bulk_reads_ch = bulk_channel
      .filter{it.technology == "single_end"}
      .map{meta -> tuple(meta,
                         file("s3://${meta.s3_prefix}/*_R1_*.fastq.gz"))}

    paired_bulk_reads_ch = bulk_channel
      .filter{it.technology == "paired_end"}
      .map{meta -> tuple(meta,
                         file("s3://${meta.s3_prefix}/*_R1_*.fastq.gz"),
                         file("s3://${meta.s3_prefix}/*_R2_*.fastq.gz"))}
    
    single_bulk
}
