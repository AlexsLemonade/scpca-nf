#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { build_celltype_ref } from './build-celltype-ref.nf'

// generate fasta and annotation files with spliced cDNA + intronic reads
process generate_reference{
  container params.SCPCATOOLS_CONTAINER
  // publish fasta and annotation files within reference directory
  publishDir params.ref_dir, mode: 'copy'
  label 'mem_32'
  maxRetries 1
  input:
    path fasta
    path gtf
    val assembly
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
  container params.SALMON_CONTAINER
  publishDir "${params.ref_dir}/salmon_index", mode: 'copy'
  label 'cpus_8'
  label 'mem_16'
  input:
    tuple path(splici_fasta), path(spliced_cdna_fasta)
    path genome
  output:
    path splici_index_dir
    path spliced_cdna_index_dir
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
  label 'cpus_12'
  label 'mem_24'
  input:
    path fasta
    path gtf
    val assembly
  output:
    path cellranger_index
  script:
    cellranger_index = "${assembly}_cellranger_full"
    """
    gunzip -c ${fasta} > genome.fasta
    gunzip -c ${gtf} > genome.gtf

    cellranger mkref \
      --genome=${cellranger_index} \
      --fasta=genome.fasta \
      --genes=genome.gtf \
      --nthreads=${task.cpus}
    """
}

process index_star{
  container params.STAR_CONTAINER
  publishDir "${params.ref_dir}/star_index", mode: 'copy'
  label 'cpus_12'
  memory '64.GB'
  input:
    path fasta
    path gtf
    val assembly
  output:
    path output_dir
  script:
    output_dir = "${assembly}.star_idx"
    """
    mkdir ${output_dir}

    # star needs uncompressed fasta & gtf
    gunzip -c ${fasta} > ${assembly}.fa
    gunzip -c ${gtf} > ${assembly}.gtf
    STAR --runMode genomeGenerate \
      --runThreadN ${task.cpus} \
      --genomeDir ${output_dir} \
      --genomeFastaFiles ${assembly}.fa \
      --genomeSAsparseD 2 \
      --sjdbGTFfile ${assembly}.gtf \
      --sjdbOverhang 100 \
      --limitGenomeGenerateRAM 64000000000

    # clean up
    rm ${assembly}.fa
    rm ${assembly}.gtf
    """
}

workflow {
  // generate splici and spliced cDNA reference fasta
  generate_reference(params.ref_fasta, params.ref_gtf, params.assembly)
  // create index using reference fastas
  salmon_index(generate_reference.out.fasta_files, params.ref_fasta)
  // create cellranger index
  cellranger_index(params.ref_fasta, params.ref_gtf, params.assembly)
  // create star index
  index_star(params.ref_fasta, params.ref_gtf, params.assembly)
}
