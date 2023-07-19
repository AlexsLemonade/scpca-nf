#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.t2g_3col_path = "s3://scpca-references/homo_sapiens/ensembl-104/annotation/Homo_sapiens.GRCh38.104.spliced_intron.tx2gene_3col.tsv"

process train_singler_models {
  container params.SCPCATOOLS_CONTAINER
  publishDir "${params.celltype_model_dir}"
  label 'cpus_4'
  label 'mem_16'
  input:
    tuple path(celltype_ref), val(ref_name)
    path tx2gene
  output:
    path celltype_model
  script:
    celltype_model = "${ref_name}_model.rds"
    """
    train_SingleR.R \
      --ref_file ${celltype_ref} \
      --output_file ${celltype_model} \
      --fry_tx2gene ${tx2gene} \
      --label_name ${params.label_name} \
      --seed ${params.seed} \
      --threads ${task.cpus}
    """
}

workflow build_celltype_ref {

  // create channel of cell type ref files and names
  celltype_refs_ch = Channel.fromPath(params.celltype_refs_metafile)
    .splitCsv(header: true, sep: '\t')
    .map{[
      celltype_ref_file = "${params.celltype_ref_dir}/${it.celltype_singler_file}",
      ref_name = it.celltype_ref_name
      ]}
    .unique() // remove any duplicates

  // train cell type references using SingleR
  train_singler_models(celltype_refs_ch, params.t2g_3col_path)
}

workflow {
  build_celltype_ref()
}
