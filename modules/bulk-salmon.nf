#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process fastp{
    container params.FASTP_CONTAINER
    label 'cpus_8'
    tag "${meta.library_id}-bulk"
    publishDir "${params.outdir}/internal/fastp"
    input: 
        tuple val(meta), path(read1), path(read2)
    output: 
        tuple val(meta), path(trimmed_reads), emit: salmon_input
        path(fastp_report), emit: report
    script: 
        trimmed_reads = "${meta.library_id}"
        fastp_report = "${meta.library_id}_fastp.html"
        """
        mkdir -p ${meta.library_id}
        cat ${read1} > ${meta.library_id}_R1_merged.fastq.gz
        cat ${read2} > ${meta.library_id}_R2_merged.fastq.gz

        fastp --in1 ${meta.library_id}_R1_merged.fastq.gz \
        ${meta.technology == 'paired_end' ? "--in2 ${meta.library_id}_R2_merged.fastq.gz" : ""} \
        --out1 ${trimmed_reads}/${meta.library_id}_R1_trimmed.fastq.gz \
        ${meta.technology == 'paired_end' ? "--out2 ${trimmed_reads}/${meta.library_id}_R2_trimmed.fastq.gz" : ""} \
        --length_required 20 \
        --thread ${task.cpus} \
        --html ${fastp_report} \
        --report_title ${meta.library_id}
        """

}

process salmon{
    container params.SALMON_CONTAINER
    label 'cpus_8'
    tag "${meta.library_id}-bulk"
    publishDir "${params.outdir}/internal/salmon/"
    input: 
        tuple val(meta), path(read_dir)
        path (index)
    output: 
        tuple val(meta), path(salmon_results)
    script:
        salmon_results = "${meta.library_id}"
        """
        salmon quant -i ${index} \
        ${meta.technology == 'paired_end' ? "-l A" : "-l U" } \
        -1 ${read_dir}/*_R1_*.fastq.gz \
        ${meta.technology == 'paired_end' ? "-2 ${read_dir}/*_R2_*.fastq.gz" : "" } \
        -o ${salmon_results} \
        --validateMappings \
        --rangeFactorizationBins 4 \
        --gcBias \
        --seqBias \
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

        fastp(bulk_reads_ch)
        salmon(fastp.out.salmon_input, params.bulk_index)
    
        emit: salmon.out
}
