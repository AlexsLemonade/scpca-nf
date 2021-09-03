#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// basic parameters
params.ref_dir = 's3://nextflow-ccdl-data/reference/homo_sapiens/ensembl-104'
params.gtf = 's3://nextflow-ccdl-data/reference/homo_sapiens/ensembl-104/annotation/Homo_sapiens.GRCh38.104.gtf.gz'
params.fasta = 's3://nextflow-ccdl-data/reference/homo_sapiens/ensembl-104/fasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz'
params.assembly = 'Homo_sapiens.GRCh38.104'

// generate fastq files with spliced cDNA + intronic reads 
process generate_fastq{
  container params.SCPCATOOLS_CONTAINTER
  // publish fasta and annotation files within reference directory 
  publishDir params.ref_dir
  input:
    path(gtf)
    path(fasta)
  output: 
    path(splici_fasta)
  script:
    splici_fasta="${params.ref_dir}/fasta/${params.assembly}.spliced_intron.txome.fa.gz"
    """
    make_splici_fasta.R \
      --gtf ${gtf} \
      --genome ${fasta}
    
    gzip annotation/*.gtf
    """
}


process salmon_index{
  container 'quay.io/biocontainers/salmon:1.4.0--hf69c8f4_0'
  publishDir "${params.ref_dir}/salmon_index", mode: 'copy'
  memory { 28.GB * task.attempt}
  cpus 8
  errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
  maxRetries 1
  input:
    path(splici_fasta)
  output:
    path(params.index_path)
  script:
    """
    salmon index \
      -t ${splici_fasta} \
      -i ${params.assembly}.spliced_intron.txome \
      -k 31 \
      -p ${task.cpus} \
    """
}


workflow {
  // generate splici reference fasta
  generate_fastq(params.gtf, params.fasta)
  // create index using splici reference fasta
  salmon_index(generate_fastq.out)
}
