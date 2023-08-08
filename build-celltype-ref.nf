#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process save_singler_refs {
  container params.SCPCATOOLS_CONTAINER
  publishDir "${params.singler_references_dir}"
  label 'mem_8'
  input:
    tuple val(ref_name), val(ref_source), path(t2g_3col_path)
  output:
    tuple val(ref_name), path(ref_file), path(t2g_3col_path)
  script:
    ref_file = "${ref_source}-${ref_name}.rds"
    """
    save_singler_refs.R \
     --ref_name ${ref_name} \
     --ref_file ${ref_file}
    """
  stub:
    ref_file = "${ref_source}-${ref_name}.rds"
    """
    touch ${ref_file}
    """

}

process train_singler_models {
  container params.SCPCATOOLS_CONTAINER
  publishDir "${params.singler_models_dir}"
  label 'cpus_4'
  label 'mem_16'
  input:
    tuple val(ref_name), path(ref_file), path(t2g_3col_path)
  output:
    path celltype_model
  script:
    celltype_model = "${ref_name}_model.rds"
    """
    train_SingleR.R \
      --ref_file ${ref_file} \
      --output_file ${celltype_model} \
      --fry_tx2gene ${t2g_3col_path} \
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

process generate_cellassign_refs {
  container params.SCPCATOOLS_CONTAINER
  publishDir "${params.cellassign_ref_dir}"
  label 'mem_8'
  input:
    tuple val(ref_name), val(ref_source), val(organs), path(ref_gtf)
    path marker_gene_file
  output:
    path ref_file
  script:
    ref_file="${ref_source}-${ref_name}.tsv"
    """
    generate_cellassign_refs.R \
      --organs "${organs}" \
      --marker_gene_file ${marker_gene_file} \
      --gtf_file ${ref_gtf} \
      --ref_mtx_file ${ref_file}
    """
  stub:
    ref_file="${ref_source}-${ref_name}.tsv"
    """
    touch ${ref_file}
    """
}

workflow build_celltype_ref {

  // read in json file with all reference paths
  ref_paths = Utils.getMetaVal(file(params.ref_json), params.celltype_organism)

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
      ref_source: it.celltype_ref_source,
      t2g_3col_path: file("${params.ref_rootdir}/${ref_paths["t2g_3col_path"]}")
      ]}

  // download and save reference files
  save_singler_refs(singler_refs_ch)

  // train cell type references using SingleR
  train_singler_models(save_singler_refs.out)

  // cellassign refs
  cellassign_refs_ch = celltype_refs_ch.cellassign
    // create a channel with ref_name, source, organs
    .map{[
      ref_name: it.celltype_ref_name,
      ref_source: it.celltype_ref_source,
      organs: it.organs,
      ref_gtf: file("${params.ref_rootdir}/${ref_paths["ref_gtf"]}")
    ]}

  generate_cellassign_refs(cellassign_refs_ch, params.panglao_marker_genes_file)

}

workflow {
  build_celltype_ref()
}
