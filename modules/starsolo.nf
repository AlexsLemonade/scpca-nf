include { index_bam } from './samtools.nf'

process starsolo {
  container params.STAR_CONTAINER
  tag "${meta.run_id}"
  label 'cpus_8'
  label 'mem_32'
  label 'disk_big'

  input:
    tuple val(meta), path(read1), path(read2), path(barcode_file), path(star_index)
  output:
    tuple val(meta), path(output_dir), emit: starsolo_dir
    tuple val(meta), path(output_bam), emit: star_bam
  script:
    tech_flag = [
      '10xv2': '',
      '10xv2_5prime': '',
      '10xv3': '--soloUMIlen 12',
      '10xv3.1': '--soloUMIlen 12',
      '10xv3_5prime': '--soloUMIlen 12',
      '10xv4': '--soloUMIlen 12'
    ]
    features_flag = meta.seq_unit == "nucleus" ? "--soloFeatures Gene GeneFull" : "--soloFeatures Gene"
    output_dir = "${meta.run_id}_star"
    output_bam = "${meta.run_id}.sorted.bam"
    """
    mkdir -p ${output_dir}/Solo.out/Gene/raw
    STAR \
      --soloType CB_UMI_Simple \
      --genomeDir ${star_index} \
      --runThreadN ${task.cpus} \
      --readFilesIn ${read2.join(',')} ${read1.join(',')} \
      --readFilesCommand gunzip -c \
      --soloCBwhitelist ${barcode_file} \
      ${tech_flag[meta.technology]} \
      ${features_flag} \
      --soloCellFilter EmptyDrops_CR \
      --outSAMtype BAM SortedByCoordinate \
      --outSAMattributes NH HI nM AS CR UR CB UB CY UY GX GN \
      --outBAMsortingThreadN 2 \
      --limitBAMsortRAM 20000000000 \
      --runDirPerm All_RWX \
      --outFileNamePrefix ${output_dir}/

    mv ${output_dir}/Aligned.sortedByCoord.out.bam ${output_bam}
    """
  stub:
    output_dir = "${meta.run_id}_star"
    output_bam = "${meta.run_id}.sorted.bam"
    """
    mkdir -p ${output_dir}
    touch ${output_bam}
    sleep 5
    """
}


workflow starsolo_map {
  take:
    singlecell_ch
    cell_barcodes

  main:
    sc_reads_ch = singlecell_ch
      .map{meta -> tuple(
        meta,
        file("${meta.files_directory}/*_R1_*.fastq.gz"),
        file("${meta.files_directory}/*_R2_*.fastq.gz"),
        file("${params.barcode_dir}/${cell_barcodes[meta.technology]}"),
        file(meta.star_index, type: 'dir')
      )}
    starsolo(sc_reads_ch)
    index_bam(starsolo.out.star_bam)

  emit:
    bam = index_bam.out
    quant = starsolo.out.starsolo_dir
}
