
process index_bam{
  container params.SAMTOOLS_CONTAINER
  input:
    tuple val(meta), path(bamfile)
  output:
    tuple val(meta), path(bamfile), path(bamfile_index)
  script:
    bamfile_index = "${bamfile}.bai"
    """
    samtools index ${bamfile} ${bamfile_index}
    """
  stub:
    bamfile_index = "${bamfile}.bai"
    """
    touch ${bamfile_index}
    """
}
