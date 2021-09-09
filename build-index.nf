#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// generate fasta and annotation files with spliced cDNA + intronic reads 
process generate_splici{
  container params.SCPCATOOLS_CONTAINER
  // publish fasta and annotation files within reference directory 
  publishDir params.ref_dir
  memory { 28.GB * task.attempt}
  errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
  maxRetries 1
  input:
    path(gtf)
    path(fasta)
    val(assembly)
  output: 
    tuple path(splici_fasta), path("annotation")
  script:
    splici_fasta="fasta/${assembly}.spliced_intron.txome.fa.gz"
    """
    make_splici_fasta.R \
      --gtf ${gtf} \
      --genome ${fasta} \
      --fasta_output fasta \
      --annotation_output annotation \
      --assembly ${assembly}
    
    gzip annotation/*.gtf
    """
}


process salmon_index{
  container 'quay.io/biocontainers/salmon:1.4.0--hf69c8f4_0'
  publishDir "${params.ref_dir}/salmon_index", mode: 'copy'
  label 'cpus_8'
  input:
    path(fasta)
  output:
    path(index_dir)
  script:
    index_dir = "${fasta}".split("\\.(fasta|fa)")[0]
    """
    salmon index \
      -t ${fasta} \
      -i ${index_dir} \
      -k 31 \
      -p ${task.cpus} \
    """
}


workflow {
  // generate splici reference fasta
  generate_splici(params.gtf, params.fasta, params.assembly)
  // create index using splici reference fasta
  salmon_index(generate_splici.out)
}
