#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.t2g_3col_path = "s3://scpca-references/homo_sapiens/ensembl-104/annotation/Homo_sapiens.GRCh38.104.spliced_intron.tx2gene_3col.tsv"

process save_singler_refs {
  container params.SCPCATOOLS_CONTAINER
  publishDir "${params.celltype_ref_dir}/singler_references"
  label 'mem_8'
  input:
    tuple val(ref_name), val(ref_database)
  output:
    tuple val(ref_name), path(ref_file)
  script:
    ref_file = "${ref_database}-${ref_name}.rds"
    """
    save_singler_refs.R \
     --ref_name ${ref_name} \
     --ref_file ${ref_file}
    """
  stub:
    ref_file = "${ref_database}-${ref_name}.rds"
    """
    touch ${ref_file}
    """

}

process train_singler_models {
  container params.SCPCATOOLS_CONTAINER
  publishDir "${params.celltype_ref_dir}/singler_models"
  label 'cpus_4'
  label 'mem_16'
  input:
    tuple val(ref_name), path(ref_file)
    path tx2gene
  output:
    path celltype_model
  script:
    celltype_model = "${ref_name}_model.rds"
    """
    train_SingleR.R \
      --ref_file ${ref_file} \
      --output_file ${celltype_model} \
      --fry_tx2gene ${tx2gene} \
      --label_name ${params.label_name} \
      --seed ${params.seed} \
      --threads ${task.cpus}
    """
  stub:
    celltype_model = "${ref_name}_model.rds"
    """
    touch ${celltype_model}
    """
}

workflow build_celltype_ref {

  // create channel of cell type ref files and names
  celltype_refs_ch = Channel.fromPath(params.celltype_ref_metadata)
    .splitCsv(header: true, sep: '\t')
    .branch{
          singler: it.celltype_method == "SingleR"
          cellassign: it.celltype_method == "CellAssign"
       }

  // singler refs to download and train
  singler_refs_ch = celltype_refs_ch.singler
    .map{[
      ref_name: it.celltype_ref_name,
      ref_database: it.celltype_ref_database
      ]}

  // download and save reference files
  save_singler_refs(singler_refs_ch)

  // train cell type references using SingleR
  train_singler_models(save_singler_refs.out, params.t2g_3col_path)

  // cellassign refs
  cellassign_refs_ch = celltype_refs_ch.cellassign
    // create a channel with ref_name, database, organs
    .map{[
      ref_name: it.celltype_ref_name,
      ref_database: it.celltype_ref_database,
      organs: it.organs
    ]}
}

workflow {
  build_celltype_ref()
}
