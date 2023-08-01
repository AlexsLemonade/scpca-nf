include { index_bam } from './samtools.nf'

process bulkmap_star{
  container params.STAR_CONTAINER
  tag "${meta.run_id}"
  label 'cpus_8'
  label 'mem_24'
  input:
    tuple val(meta), path(read1), path(read2), path(star_index)
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
  stub:
    output_bam = "${meta.run_id}.sorted.bam"
    """
    touch ${output_bam}
    """
}

workflow star_bulk{
  take:
    bulk_channel // a channel with a map of metadata for each rna library to process

  main:
    // create tuple of (metadata map, [Read 1 files], [Read 2 files])
    bulk_reads_ch = bulk_channel
        .map{meta -> tuple(meta,
                           file("${meta.files_directory}/*_{R1,R1_*}.fastq.gz"),
                         file("${meta.files_directory}/*_{R2,R2_*}.fastq.gz"),
                           file("${meta.star_index}"))}
    // map and index
    bulkmap_star(bulk_reads_ch) \
    | index_bam

  emit:
    index_bam.out
}
