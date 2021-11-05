#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process fastp{
    container params.FASTP_CONTAINER
    label 'cpus_8'
    tag "${meta.library_id}-bulk"
    input: 
        tuple val(meta), path(read1), path(read2)
    output: 
        tuple val(meta), path(trimmed_reads)
    script: 
        trimmed_reads = "${meta.library_id}"
        fastp_report = "${meta.library_id}_fastp.html"
        """
        mkdir -p ${trimmed_reads}
        fastp --in1 <(gunzip -c ${read1}) --out1 ${trimmed_reads}/${meta.library_id}_R1_trimmed.fastq.gz \
        ${meta.technology == 'paired_end' ? "--in2 <(gunzip -c ${read2}) --out2 ${trimmed_reads}/${meta.library_id}_R2_trimmed.fastq.gz" : ""} \
        --length_required 20 \
        --thread ${task.cpus}
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
        -l A \
        ${meta.technology == 'paired_end' ? "-1": "-r"} ${read_dir}/*_R1_*.fastq.gz \
        ${meta.technology == 'paired_end' ? "-2 ${read_dir}/*_R2_*.fastq.gz" : "" } \
        -o ${salmon_results} \
        --validateMappings \
        --rangeFactorizationBins 4 \
        --gcBias \
        --seqBias \
        --threads ${task.cpus}
        """

}

process group_tximport {
    container params.SCPCATOOLS_CONTAINER
    publishDir "${params.outdir}/publish/${project_id}"
    input:
        tuple val(project_id), path(salmon_directories)
        path(tx2gene)
    output:
        path(tximport_file)
    script:
        tximport_file = "${project_id}_bulk_quant.tsv"
        """
        ls -d ${salmon_directories} > salmon_directories.txt

        merge_counts_tximport.R \
          --project_id ${project_id} \
          --salmon_dir salmon_directories.txt \
          --output_file ${tximport_file} \
          --tx2gene ${tx2gene}
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
        salmon(fastp.out, params.bulk_index)

        grouped_salmon_ch = salmon.out
            .map{[it[0]["project_id"], it[1]]}
            .groupTuple(by: 0)

        group_tximport(grouped_salmon_ch, params.t2g_bulk_path)
    emit: 
        group_tximport.out
}
