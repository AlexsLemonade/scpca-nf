include { index_bam } from './samtools.nf'

process bulkmap_star {
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

workflow star_bulk {
  take:
    bulk_channel // a channel with a map of metadata for each rna library to process

  main:
    // create list of (metadata map, [Read 1 files], [Read 2 files], star index)
    // regex to ensure correct file names if R1 or R2 are in sample identifier
    bulk_reads_ch = bulk_channel
      .map{ meta ->
        def fastq_files = files("${meta.files_directory}/*.fastq.gz", checkIfExists: true)
        // add R1 and R2 regex to ensure correct file names if R1 or R2 are in sample identifier
        def R1_files = fastq_files.findAll{ it.name =~ /_R1(_\d+)?\.fastq\.gz$/ }
        def R2_files = fastq_files.findAll{ it.name =~ /_R2(_\d+)?\.fastq\.gz$/ }

        // check that appropriate files were found
        assert R1_files: "No R1 files were found in ${meta.files_directory}."
        if (meta.technology == 'paired_end') {
          assert R2_files: "No R2 files were found in ${meta.files_directory}."
        }

        def star_index = file(meta.star_index, type: 'dir', checkIfExists: true)
        [meta, R1_files, R2_files, star_index]
      }

    // map and index
    bulkmap_star(bulk_reads_ch) \
      | index_bam

  emit:
    index_bam.out
}
