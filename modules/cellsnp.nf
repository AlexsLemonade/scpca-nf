
process cellsnp{
  container params.CELLSNP_CONTAINER
  label 'cpus_8'
  label 'mem_16'
  tag "${meta_star.run_id}"
  input:
    tuple val(meta_star), path(star_bam), path(star_bai), path(star_quant),
          val(meta_mpileup), path(vcf_file)
  output:
    tuple val(meta), path(outdir), path(vcf_file)
  script:
    meta = meta_star
    meta.sample_ids = meta_mpileup.sample_ids
    meta.bulk_run_ids = meta_mpileup.bulk_run_ids
    quant_dir = meta_star.seq_unit == "nucleus" ? "GeneFull" : "Gene"
    barcodes = "${star_quant}/Solo.out/${quant_dir}/filtered/barcodes.tsv"
    outdir = "${meta.run_id}-cellSNP"
    """
    cellsnp-lite \
      --samFile ${star_bam} \
      --barcodeFile ${barcodes} \
      --regionsVCF <(gunzip -c ${vcf_file}) \
      --nproc ${task.cpus} \
      --outDir ${outdir} \
      --minMAF=0.1 \
      --minCOUNT=20 \
      --gzip
    """
}

process vireo{
  container params.VIREO_CONTAINER
  publishDir "${meta.vireo_publish_dir}", mode: 'copy'
  tag "${meta.run_id}"
  label 'cpus_8'
  label 'mem_16'
  input:
    tuple val(meta), path(cellsnp_dir), path(vcf_file)
  output:
    tuple val(meta), path(outdir)
  script:
    outdir = file(meta.vireo_dir).name
    meta_json = Utils.makeJson(meta)
    """
    vireo \
      --cellData ${cellsnp_dir} \
      --donorFile ${vcf_file}  \
      --outDir ${outdir} \
      --nproc ${task.cpus}

    echo '${meta_json}' > ${outdir}/scpca-meta.json
    """
}

workflow cellsnp_vireo {
  take:
    starsolo_bam_ch //channel of [meta, bamfile, bam.bai]
    starsolo_quant_ch //channel of [meta, starsolo_dir]
    mpileup_vcf_ch // channel of [meta_mpileup, vcf_file]
  main:
    mpileup_ch = mpileup_vcf_ch
      .map{[it[0].multiplex_library_id] + it} // pull out library id for combining
    star_mpileup_ch = starsolo_bam_ch.map{[it[0].library_id] + it} // add library id at start
      .join(starsolo_quant_ch.map{[it[0].library_id] + it}, by: 0, failOnDuplicate: true, failOnMismatch: true) // join starsolo outs by library_id
      .map{[it[0], it[1], it[2], it[3], it[5]]} // remove redundant meta
      .join(mpileup_ch, by: 0, failOnDuplicate: true, failOnMismatch: true) // join starsolo and mpileup by library id
      .map{it.drop(1)} // drop library id
      //result: [meta, star_bam, star_bai, star_quant, meta_mpileup, vcf_file]

    cellsnp(star_mpileup_ch) \
      | vireo
  emit:
    vireo.out
}

