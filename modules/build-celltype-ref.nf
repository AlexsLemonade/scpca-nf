#!/usr/bin/env nextflow
nextflow.enable.dsl=2

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
      --seed ${params.seed} \
      --threads ${task.cpus}
    """
}

workflow build_celltype_ref {
  take: celltype_refs_file

  main:
    // create channel of cell type ref files and names
    celltype_refs_ch = celltype_refs_file
        .map{[
        celltype_ref_file = "${params.celltype_ref_dir}/${it.filename}",
        ref_name = it.reference
        ]}

    // train cell type references using SingleR
    train_singler_models(celltype_refs_ch, params.t2g_3col_path)
}
