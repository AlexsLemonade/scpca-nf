#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// generate fasta and annotation files with spliced cDNA + intronic reads 
process generate_reference{
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
    tuple path(splici_fasta), path(spliced_cdna_fasta), emit: fasta_files
    tuple path("annotation/*.gtf.gz"), path("annotation/*.tsv"), path("annotation/*.txt"),  emit: annotations
  script:
    splici_fasta="fasta/${assembly}.spliced_intron.txome.fa.gz"
    spliced_cdna_fasta="fasta/${assembly}.spliced_cdna.txome.fa.gz"
    """
    make_reference_fasta.R \
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
    tuple path(splici_fasta), path(spliced_cdna_fasta)
    path(genome)
  output:
    path(splici_index_dir)
    path(spliced_cdna_index_dir)
  script:
    splici_index_dir = "${splici_fasta}".split("\\.(fasta|fa)")[0]
    spliced_cdna_index_dir = "${spliced_cdna_fasta}".split("\\.(fasta|fa)")[0]
    """
    salmon index \
      -t ${splici_fasta} \
      -i ${splici_index_dir} \
      -k 31 \
      -p ${task.cpus} \

    gunzip -c ${genome} \
      |grep "^>" | cut -d " " -f 1 \
      |sed -e 's/>//g' > decoys.txt
    cat ${spliced_cdna_fasta} ${genome} > gentrome.fa.gz

    salmon index \
      -t gentrome.fa.gz \
      -d decoys.txt \
      -i ${spliced_cdna_index_dir} \
      -k 31 \
      -p ${task.cpus} \
    """
}

process cellranger_index{
  container params.CELLRANGER_CONTAINER
  publishDir "${params.ref_dir}/cellranger_index", mode: 'copy'
  label 'cpus_8'
  input:
    path(fasta)
    path(gtf)
    val(assembly)
  output:
    path(cellranger_index)
  script:
    cellranger_index = "${assembly}_cdna_cellranger"
    """
    gunzip -c ${fasta} > genome.fasta
    gunzip -c ${gtf} > genome.gtf
    cellranger mkgtf \
      genome.gtf \
      filtered.gtf \
      --attribute=gene_biotype:protein_coding
    
    cellranger mkref \
      --genome=${cellranger_index} \
      --fasta=genome.fasta \
      --genes=filtered.gtf \
      --nthreads=${task.cpus}
    """
}



workflow {
  // generate splici and spliced cDNA reference fasta
  generate_reference(params.gtf, params.fasta, params.assembly)
  // create index using reference fastas
  salmon_index(generate_reference.out.fasta_files, params.fasta)
  // create cellranger index 
  cellranger_index(params.fasta, params.gtf, params.assembly)
}
