#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// parameters
params.assembly  = 'Homo_sapiens.GRCh38.104'
params.ref_dir   = 's3://nextflow-ccdl-data/reference/homo_sapiens/ensembl-104'
params.ref_fasta = 's3://nextflow-ccdl-data/reference/homo_sapiens/ensembl-104/fasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz'
params.ref_gtf   = 's3://nextflow-ccdl-data/reference/homo_sapiens/ensembl-104/annotation/Homo_sapiens.GRCh38.104.gtf.gz'

STARCONTAINER = 'quay.io/biocontainers/star:2.7.9a--h9ee0642_0'

process index_star{
  container STARCONTAINER
  publishDir "${params.ref_dir}/star_index", mode: 'copy'
  memory "64.GB"
  cpus "8"
  input:
    path ref_fasta
    path ref_gtf
  output:
    path output_dir
  script:
    output_dir = "${params.assembly}.star_idx"
    """
    mkdir ${output_dir}

    # star needs uncompressed fasta & gtf
    gunzip -c ${ref_fasta} > ${params.assembly}.fa
    gunzip -c ${ref_gtf} > ${params.assembly}.gtf
    STAR --runMode genomeGenerate \
      --runThreadN ${task.cpus} \
      --genomeDir ${output_dir} \
      --genomeFastaFiles ${params.assembly}.fa \
      --sjdbGTFfile ${params.assembly}.gtf \
      --sjdbOverhang 100 \
      --limitGenomeGenerateRAM 64000000000
    
    # clean up
    rm ${params.assembly}.fa
    rm ${params.assembly}.gtf
    """
}

workflow{
  // index with star
  index_star(params.ref_fasta, params.ref_gtf)
}
