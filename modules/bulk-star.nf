#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process bulkmap_star{
  container params.STAR_CONTAINER
  tag "${meta.run_id}"
  memory "32.GB"
  cpus "8"
  input:
    tuple val(meta), path(read1), path(read2)
    path star_index
  output:
    tuple val(meta), path(output_bam)
  script:
    output_bam = "${meta.run_id}.sorted.bam"
    """
    STAR \
      --genomeDir ${star_index} \
      --runThreadN ${task.cpus} \
      --readFilesIn ${read1.join(',')} \
      ${meta.technology == 'paired_end' ? read2.join(',') : ""} \
      --readFilesCommand gunzip -c \
      --outSAMtype BAM SortedByCoordinate

    mv Aligned.sortedByCoord.out.bam ${output_bam}
    """
}

process index_bam{
  container params.SAMTOOLS_CONTAINER
  tag "${meta.run_id}"
  input:
    tuple val(meta), path(bamfile)
  output:
    tuple val(meta), path(bamfile), path(bamfile_index)
  script:
    bamfile_index = "${bamfile}.bai"
    """
    samtools index ${bamfile} ${bamfile_index}
    """
}

workflow star_bulk{
  take: 
    bulk_channel // a channel with a map of metadata for each rna library to process
    
  main: 
    // create tuple of (metadata map, [Read 1 files], [Read 2 files])
    bulk_reads_ch = bulk_channel
        .map{meta -> tuple(meta,
                            file("s3://${meta.s3_prefix}/*_R1_*.fastq.gz"),
                            file("s3://${meta.s3_prefix}/*_R2_*.fastq.gz"))}
    // map and index
    bulkmap_star(bulk_reads_ch, params.star_index) \
    | index_bam
  
  emit: 
    index_bam.out
}
