#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process fastp{
    container params.FASTP_CONTAINER
    label 'cpus_8'
    tag "${meta.library_id}-bulk"
    publishDir "${params.outdir}/internal/fastp"
    input: 
        tuple val(meta),path(read1), path(read2)
    output: 
        tuple val(meta),path(trimmed_reads)
    script: 
        trimmed_reads = "${meta.library_id}"
        """
        mkdir -p ${meta.library_id}

        fastp --in1 ${read1} \
        ${meta.technology == 'paired_end' ? "--in2 ${read2}":""} \
        --out1 ${trimmed_reads}/${meta.library_id}-trimmed-R1.fastq.gz \
        ${meta.technology == 'paired_end' ? "--out2 ${trimmed_reads}/${meta.library_id}-trimmed-R2.fastq.gz":""} \
        --length_required 20 \
        --html ${meta.library_id}_fastp.html \
        --json ${meta.library_id}_fastp.json \
        --report_title ${meta.library_id}
        """

}

process salmon{
    container params.SALMON_CONTAINER
    label 'cpus_8'
    tag "${meta.library_id}-bulk"
    publishDir "${params.outdir}/internal/salmon/${meta.library_id}"
    input: 
        tuple val(meta),path(trimmed_reads)
    output: 
        tuple val(meta),path(salmon_results)
    script:
        salmon_results = "${meta.library_id}-salmon"
        """
        salmon quant -i ${params.bulk_index} /
        -libType A /
        -1 ${trimmed_reads}/${meta.library_id}-trimmed-R1.fastq.gz /
        ${meta.technology == 'paired_end' ? "-2 ${trimmed_reads}/${meta.library_id}-trimmed-R2.fastq.gz":""} /
        -o ${salmon_results} /
        --validateMappings /
        --rangeFactorizationBins 4 /
        --gcBias /
        --seqBias /
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
