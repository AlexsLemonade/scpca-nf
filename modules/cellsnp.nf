
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
  container params.CONDA_CONTAINER
  publishDir "${meta.vireo_publish_dir}"
  tag "${meta.run_id}"
  label 'cpus_8'
  label 'mem_16'
  input:
    tuple val(meta), path(cellsnp_dir), path(vcf_file)
  output:
    tuple val(meta), path(outdir)
  script:
    outdir = file(meta.vireo_dir).name
    """
    pip install vireoSNP==0.5.6
    vireo \
      --cellData ${cellsnp_dir} \
      --donorFile ${vcf_file}  \
      --outDir ${outdir} \
      --nproc ${task.cpus}
    """
    // write out meta object
    meta_file = file("${outdir}/scpca.json")
    write_meta(meta, meta_file)
}


