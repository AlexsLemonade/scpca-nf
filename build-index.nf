#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { build_celltype_ref } from './build-celltype-ref.nf'

// generate fasta and annotation files with spliced cDNA + intronic reads
process generate_reference {
  container params.SCPCATOOLS_CONTAINER
  // publish fasta and annotation files within reference directory
  publishDir "${params.ref_rootdir}/${meta.ref_dir}", mode: 'copy'
  label 'mem_32'
  maxRetries 1
  input:
    tuple val(ref_name), val(meta), path(fasta), path(gtf)
  output:
    tuple val(ref_name), val(meta), emit: ref_info
    tuple path(splici_fasta), path(spliced_cdna_fasta), emit: fasta_files
    tuple path("annotation/*.gtf.gz"), path("annotation/*.tsv"), path("annotation/*.txt"),  emit: annotations
  script:
    splici_fasta = "fasta/" + file("${meta.splici_index}").name + ".fa.gz"
    spliced_cdna_fasta =  "fasta/" + file("${meta.salmon_bulk_index}").name + ".fa.gz"
    """
    make_reference_fasta.R \
      --gtf ${gtf} \
      --genome ${fasta} \
      --fasta_output fasta \
      --annotation_output annotation \
      --reference_name ${ref_name}

    gzip annotation/*.gtf
    """
}


process salmon_index {
  container params.SALMON_CONTAINER
  publishDir "${params.ref_rootdir}/${meta.ref_dir}/salmon_index", mode: 'copy'
  label 'cpus_8'
  label 'mem_16'
  input:
    tuple path(splici_fasta), path(spliced_cdna_fasta)
    tuple val(ref_name), val(meta), path(fasta), path(gtf)
  output:
    path splici_index_dir
    path spliced_cdna_index_dir
  script:
    splici_index_dir = file("${meta.splici_index}").name
    spliced_cdna_index_dir = file("${meta.salmon_bulk_index}").name
    """
    salmon index \
      -t ${splici_fasta} \
      -i ${splici_index_dir} \
      -k 31 \
      -p ${task.cpus} \

    gunzip -c ${fasta} \
      | grep "^>" | cut -d " " -f 1 \
      | sed -e 's/>//g' > decoys.txt
    cat ${spliced_cdna_fasta} ${fasta} > gentrome.fa.gz

    salmon index \
      -t gentrome.fa.gz \
      -d decoys.txt \
      -i ${spliced_cdna_index_dir} \
      -k 31 \
      -p ${task.cpus} \
    """
}

process cellranger_index {
  container params.CELLRANGER_CONTAINER
  publishDir "${params.ref_rootdir}/${meta.ref_dir}/cellranger_index", mode: 'copy'
  label 'cpus_12'
  label 'mem_24'
  input:
    tuple val(ref_name), val(meta), path(fasta), path(gtf)
  output:
    path cellranger_index
  script:
    cellranger_index = file("${meta.cellranger_index}").name
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

process star_index {
  container params.STAR_CONTAINER
  publishDir "${params.ref_rootdir}/${meta.ref_dir}/star_index", mode: 'copy'
  label 'cpus_12'
  memory '64.GB'
  input:
    tuple val(ref_name), val(meta), path(fasta), path(gtf)
  output:
    path output_dir
  script:
    output_dir = file("${meta.star_index}").name
    """
    mkdir ${output_dir}

    # star needs uncompressed fasta & gtf
    gunzip -c ${fasta} > ${ref_name}.fa
    gunzip -c ${gtf} > ${ref_name}.gtf
    STAR --runMode genomeGenerate \
      --runThreadN ${task.cpus} \
      --genomeDir ${output_dir} \
      --genomeFastaFiles ${ref_name}.fa \
      --genomeSAsparseD 2 \
      --sjdbGTFfile ${ref_name}.gtf \
      --sjdbOverhang 100 \
      --limitGenomeGenerateRAM 64000000000

    # clean up
    rm ${ref_name}.fa
    rm ${ref_name}.gtf
    """
}

workflow {

  // read in json file with all reference paths
  ref_paths = Utils.readMeta(file(params.ref_json))

  // read in metadata with all organisms to create references for
  ref_ch = Channel.fromPath(params.ref_metadata)
    .splitCsv(header: true, sep: '\t')
    .map{
      def reference_name = "${it.organism}.${it.assembly}.${it.version}";
      // reference name & reference file paths for each organism
      [reference_name, ref_paths[reference_name]]
    }
    .map{it +[
      file("${params.ref_rootdir}/${it[1]["ref_fasta"]}"), // path to fasta
      file("${params.ref_rootdir}/${it[1]["ref_gtf"]}") // path to gtf
    ]}


  // generate splici and spliced cDNA reference fasta
  generate_reference(ref_ch)
  // create index using reference fastas
  salmon_index(generate_reference.out.fasta_files, ref_ch)
  // create cellranger index
  cellranger_index(ref_ch)
  // create star index
  star_index(ref_ch)

}
