#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include {  getMetaVal; pullthroughContainer } from './lib/utils.nf'

process save_singler_refs {
  container pullthroughContainer(params.scpcatools_container, params.pullthrough_registry)
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
    touch "${ref_name}_${ref_source}_0-0-0.rds"
    """

}

process train_singler_models_transcriptome {
  container pullthroughContainer(params.scpcatools_container, params.pullthrough_registry)
  publishDir "${params.singler_models_dir}"
  label 'cpus_4'
  label 'mem_16'
  input:
    tuple val(ref_name), path(ref_file)
    path t2g_3col_path
    val ref_assembly // corresponds to assembly in meta.json files output from main.nf
  output:
    path celltype_model
  script:
    gene_set_version = ref_assembly.tokenize('.')
      .takeRight(2) // take the last two elements which have assembly and version
      .join('-') // join to get GRCh38-104
    date_str = java.time.LocalDate.now().toString() // get current date in ISO8601 format
    celltype_model = "${ref_file.baseName}_${gene_set_version}_${date_str}_model.rds"
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
    gene_set_version = ref_assembly.tokenize('.')
      .takeRight(2)
      .join('-')
    date_str = java.time.LocalDate.now().toString()
    celltype_model = "${ref_file.baseName}_${gene_set_version}_${date_str}_model.rds"
    """
    touch ${celltype_model}
    """
}

process train_singler_models_flex {
  container pullthroughContainer(params.scpcatools_container, params.pullthrough_registry)
  publishDir "${params.singler_models_dir}"
  label 'cpus_4'
  label 'mem_16'
  input:
    tuple val(ref_name), path(ref_file)
    path flex_probeset
    val probeset_version
  output:
    path celltype_model
  script:
    date_str = java.time.LocalDate.now().toString() // get current date in ISO8601 format
    celltype_model = "${ref_file.baseName}_${probeset_version}_${date_str}_model.rds"
    """
    train_SingleR.R \
      --ref_file ${ref_file} \
      --output_file ${celltype_model} \
      --flex_probeset ${flex_probeset} \
      --label_name ${params.singler_label_name} \
      --seed ${params.seed} \
      --threads ${task.cpus}
    """
  stub:
    date_str = java.time.LocalDate.now().toString()
    celltype_model = "${ref_file.baseName}_${probeset_version}_${date_str}_model.rds"
    """
    touch ${celltype_model}
    """
}

process catalog_singler_models {
  container pullthroughContainer(params.tidyverse_container, params.pullthrough_registry)
  publishDir "${params.singler_models_dir}"
  input:
    val celltype_references
  output:
    path "singler_models.tsv"
  script:
    """
    make_celltype_ref_table.R "${celltype_references}" singler_models.tsv
    """
  stub:
    """
    touch singler_models.tsv
    """
}

process generate_cellassign_refs {
  container pullthroughContainer(params.scpcatools_container, params.pullthrough_registry)
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

process catalog_cellassign_refs {
  container pullthroughContainer(params.tidyverse_container, params.pullthrough_registry)
  publishDir "${params.cellassign_ref_dir}"
  input:
    val celltype_references
  output:
    path "cellassign_references.tsv"
  script:
    """
    make_celltype_ref_table.R "${celltype_references}" cellassign_references.tsv
    """
  stub:
    """
    touch cellassign_references.tsv
    """
}

workflow build_celltype_ref {

  // read in json file with all reference paths
  ref_paths = getMetaVal(file(params.ref_json), params.celltype_organism)
  // get path to tx2gene and gtf
  t2g_3col_path = file("${params.ref_rootdir}/${ref_paths["t2g_3col_path"]}")
  ref_gtf = file("${params.ref_rootdir}/${ref_paths["ref_gtf"]}")

  // create channel of cell type ref files and names
  celltype_refs_ch = channel.fromPath(params.celltype_ref_metadata)
    .splitCsv(header: true, sep: '\t')
    .branch{ it ->
      singler: it.celltype_method == "SingleR"
      cellassign: it.celltype_method == "CellAssign"
    }


  // singler refs to download and train
  singler_refs_ch = celltype_refs_ch.singler
    .map{ ref ->
      [ref.celltype_ref_name, ref.celltype_ref_source]
    }

  // download and save reference files
  save_singler_refs(singler_refs_ch)

  // train cell type references using SingleR with transcriptome
  train_singler_models_transcriptome(save_singler_refs.out, t2g_3col_path, params.celltype_organism)

  // train cell type references with probe sets for 10x flex
  train_singler_models_flex(save_singler_refs.out, file(params.flex_probeset_file), params.celltype_probeset_version)

  // combine all output model files and join join file names into a comma separated string
  singler_models = train_singler_models_transcriptome.out.mix(train_singler_models_flex.out)
    .reduce{ a, b -> "$a,$b" }
  catalog_singler_models(singler_models)

  // cellassign refs
  cellassign_refs_ch = celltype_refs_ch.cellassign
    // create a channel with ref_name, source, organs
    .map{ ref ->
      [ref.ref_name, ref.ref_source, ref.organs]
    }

  generate_cellassign_refs(cellassign_refs_ch, ref_gtf, params.panglao_marker_genes_file)

  // join reference file names into a comma separated string
  cellassign_refs = generate_cellassign_refs.out.reduce{a, b -> "$a,$b"}
  catalog_cellassign_refs(cellassign_refs)
}

workflow {
  build_celltype_ref()
}
