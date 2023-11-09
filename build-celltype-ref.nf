#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process save_singler_refs {
  container params.SCPCATOOLS_CONTAINER
  publishDir "${params.singler_references_dir}"
  label 'mem_8'
  input:
    tuple val(ref_name), val(ref_source)
  output:
    tuple val(ref_name), path("${ref_name}_${ref_source}_*.rds")
  script:
    """
    save_singler_refs.R \
     --ref_name ${ref_name} \
     --ref_file_prefix "${ref_name}_${ref_source}"
    """
  stub:
    // fill in a dummy version since we grab that as part of the script
    """
    touch "${ref_name}_${ref_source}_v0-0-0.rds"
    """

}

process train_singler_models {
  container params.SCPCATOOLS_CONTAINER
  publishDir "${params.singler_models_dir}"
  label 'cpus_4'
  label 'mem_16'
  input:
    tuple val(ref_name), path(ref_file)
    path t2g_3col_path
  output:
    path celltype_model
  script:
    ref_file_basename = file("${ref_file}").baseName
    celltype_model = "${ref_file_basename}_model.rds"
    """
    train_SingleR.R \
      --ref_file ${ref_file} \
      --output_file ${celltype_model} \
      --fry_tx2gene ${t2g_3col_path} \
      --label_name ${params.singler_label_name} \
      --seed ${params.seed} \
      --threads ${task.cpus}
    """
  stub:
    ref_file_basename = file("${ref_file}").baseName
    celltype_model = "${ref_file_basename}_model.rds"
    """
    touch ${celltype_model}
    """
}

process generate_cellassign_refs {
  container params.SCPCATOOLS_CONTAINER
  publishDir "${params.cellassign_ref_dir}"
  label 'mem_8'
  input:
    tuple val(ref_name), val(ref_source), val(organs)
    path ref_gtf
    path marker_gene_file
  output:
    path ref_file
  script:
    // get ref version from filename
    // this requires the date stored in the filename to be in ISO8601 format
    ref_version = (marker_gene_file =~ /.+(20[0-9]{2}\-[0-9]{2}\-[0-9]{2}).tsv/)[0][1]
    ref_file = "${ref_name}_${ref_source}_${ref_version}.tsv"
    """
    generate_cellassign_refs.R \
      --organs "${organs}" \
      --marker_gene_file ${marker_gene_file} \
      --gtf_file ${ref_gtf} \
      --ref_mtx_file ${ref_file}
    """
  stub:
    ref_version = (marker_gene_file =~ /.+(20[0-9]{2}\-[0-9]{2}\-[0-9]{2}).tsv/)[0][1]
    ref_file = "${ref_name}_${ref_source}_${ref_version}.tsv"
    """
    touch ${ref_file}
    """
}

workflow build_celltype_ref {

  // read in json file with all reference paths
  ref_paths = Utils.getMetaVal(file(params.ref_json), params.celltype_organism)
  // get path to tx2gene and gtf
  t2g_3col_path = file("${params.ref_rootdir}/${ref_paths["t2g_3col_path"]}")
  ref_gtf = file("${params.ref_rootdir}/${ref_paths["ref_gtf"]}")

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
      ref_source: it.celltype_ref_source
      ]}

  // download and save reference files
  save_singler_refs(singler_refs_ch)

  // train cell type references using SingleR
  train_singler_models(save_singler_refs.out, t2g_3col_path)

  // cellassign refs
  cellassign_refs_ch = celltype_refs_ch.cellassign
    // create a channel with ref_name, source, organs
    .map{[
      ref_name: it.celltype_ref_name,
      ref_source: it.celltype_ref_source,
      organs: it.organs
    ]}

  generate_cellassign_refs(cellassign_refs_ch, ref_gtf, params.panglao_marker_genes_file)

}

workflow {
  build_celltype_ref()
}
