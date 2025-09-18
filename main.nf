#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// include processes from modules
include { map_quant_rna } from './modules/af-rna.nf'
include { map_quant_feature } from './modules/af-features.nf'
include { bulk_quant_rna } from './modules/bulk-salmon.nf'
include { genetic_demux_vireo } from './modules/genetic-demux.nf'
include { spaceranger_quant } from './modules/spaceranger.nf'
include { flex_quant } from './modules/cellranger-flex.nf'
include { generate_sce; generate_sce_with_feature; generate_sce_cellranger; cellhash_demux_sce; genetic_demux_sce; post_process_sce} from './modules/sce-processing.nf'
include { cluster_sce } from './modules/cluster-sce.nf'
include { annotate_celltypes } from './modules/classify-celltypes.nf'
include { call_cnvs } from './modules/call-cnvs.nf'
include { qc_publish_sce } from './modules/publish-sce.nf'
include { sce_to_anndata } from './modules/export-anndata.nf'

def check_parameters() {
  // parameter check function
  def param_error = false
  if (!file(params.run_metafile).exists()) {
    log.error("The 'run_metafile' file '${params.run_metafile}' can not be found.")
    param_error = true
  }

  if (!file(params.sample_metafile).exists()) {
    log.error("The 'sample_metafile' file '${params.sample_metafile}' can not be found.")
    param_error = true
  }

  def resolution_strategies = ['cr-like', 'full', 'cr-like-em', 'parsimony', 'trivial']
  if (!resolution_strategies.contains(params.af_resolution)) {
    log.error("'af_resolution' must be one of the following: ${resolution_strategies}")
    param_error = true
  }

  if (params.cellhash_pool_file && !file(params.cellhash_pool_file).exists()){
    log.error("The 'cellhash_pool_file' file ${params.cellhash_pool_file} can not be found.")
    param_error = true
  }

  perform_celltyping = params.perform_celltyping // create separate variable to use instead
  // CNV inference checks
  if (params.perform_cnv_inference) {
    if (!perform_celltyping) {
      log.warn("To call CNVs, cell typing must be performed as well. Setting `--perform_celltyping` to true")
      perform_celltyping = true
    }
    if (!file(params.diagnosis_groups_file).exists()) {
      log.error("The 'diagnosis_groups_file' file '${params.diagnosis_groups_file}' can not be found.")
      param_error = true
    }
    if (!file(params.diagnosis_celltypes_file).exists()) {
      log.error("The 'diagnosis_celltypes_file' file '${params.diagnosis_celltypes_file}' can not be found.")
      param_error = true
    }
  }

  // cell type annotation file checks
  if (perform_celltyping) {
    if (!file(params.project_celltype_metafile).exists()) {
      log.error("The 'project_celltype_metafile' file '${params.project_celltype_metafile}' can not be found.")
      param_error = true
    }
    if (!file(params.celltype_ref_metadata).exists()) {
      log.error("The 'celltype_ref_metadata' file '${params.celltype_ref_metadata}' can not be found.")
      param_error = true
    }
    // check that reference files related to consensus cell types exist
    if (!file(params.consensus_ref_file).exists()) {
      log.error("The 'consensus_ref_file' file '${params.consensus_ref_file}' can not be found.")
      param_error = true
    }
    if (!file(params.validation_groups_file).exists()) {
      log.error("The 'validation_groups_file' file '${params.validation_groups_file}' can not be found.")
      param_error = true
    }
    if (!file(params.validation_markers_file).exists()) {
      log.error("The 'validation_markers_file' file '${params.validation_markers_file}' can not be found.")
      param_error = true
    }
    if (!file(params.validation_palette_file).exists()) {
      log.error("The 'validation_palette_file' file '${params.validation_palette_file}' can not be found.")
      param_error = true
    }
  }

  if (param_error) {
    System.exit(1)
  }
}



// Main workflow
workflow {

  // check parameters before starting workflow
  check_parameters()

  /// Define setup variables
  // 10X barcode files
  def cell_barcodes = [
    '10Xv2': '737K-august-2016.txt',
    '10Xv2_5prime': '737K-august-2016.txt',
    '10Xv3': '3M-february-2018.txt',
    '10Xv3.1': '3M-february-2018.txt',
    '10Xv3_5prime': '3M-5pgex-jan-2023.txt.gz',
    '10Xv4': '3M-3pgex-may-2023_TRU.txt.gz',
    'CITEseq_10Xv2': '737K-august-2016.txt',
    'CITEseq_10Xv3': '3M-february-2018.txt',
    'CITEseq_10Xv3.1': '3M-february-2018.txt',
    'cellhash_10Xv2': '737K-august-2016.txt',
    'cellhash_10Xv3': '3M-february-2018.txt',
    'cellhash_10Xv3.1': '3M-february-2018.txt'
  ]

  // 10X flex probe set files
  def flex_probesets = [
    '10Xflex_v1.1_single': 'Chromium_Human_Transcriptome_Probe_Set_v1.1.0_GRCh38-2024-A.csv',
    '10Xflex_v1.1_multi': 'Chromium_Human_Transcriptome_Probe_Set_v1.1.0_GRCh38-2024-A.csv'
  ]

  // supported technologies
  def single_cell_techs = cell_barcodes.keySet()
  def flex_techs = flex_probesets.keySet()
  def bulk_techs = ['single_end', 'paired_end']
  def spatial_techs = ['visium']
  def all_techs = single_cell_techs + bulk_techs + spatial_techs
  def rna_techs = single_cell_techs.findAll{it.startsWith('10Xv')}
  def citeseq_techs = single_cell_techs.findAll{it.startsWith('CITEseq')}
  def cellhash_techs = single_cell_techs.findAll{it.startsWith('cellhash')}

  // used when a given file is not defined
  def empty_file = "${projectDir}/assets/NO_FILE"

  // select runs to use
  def run_ids = []
  def project_ids = []
  if (params.project) {
    // projects will use all runs in the project & supersede run_ids
    // allow for processing of multiple projects at once
    project_ids = params.project?.tokenize(',') ?: []
  } else {
    run_ids = params.run_ids?.tokenize(',') ?: []
  }
  def run_all = run_ids[0] == "All"
  if (run_all) {
    log.info("Executing workflow for all runs in the run metafile.")
  }

  def ref_paths = Utils.readMeta(file(params.ref_json))

  all_runs_ch = Channel.fromPath(params.run_metafile)
    .splitCsv(header: true, sep: '\t')
    // only technologies we know how to process
    .filter{it.technology in all_techs}
    // use only the rows in the run_id list (run, library, or sample can match)
    // or run by project or submitter if the project parameter is set
    .filter{
      run_all
      || (it.scpca_run_id in run_ids)
      || (it.scpca_library_id in run_ids)
      || (it.scpca_sample_id in run_ids)
      || (it.submitter in project_ids)
      || (it.scpca_project_id in project_ids)
    }
    .branch{ it ->
      known_ref: it.sample_reference in ref_paths
      unknown_ref: true
    }

  // warn about unknown references for any samples
  all_runs_ch.unknown_ref.subscribe{ it ->
    log.warn("The sample reference '${it.sample_reference}' for run '${it.scpca_run_id}' is not known and this run will not be processed.")
  }

  // convert row data to a metadata map, keeping columns we will need (& some renaming) and reference paths
  unfiltered_runs_ch = all_runs_ch.known_ref
    .map{
      def sample_refs = ref_paths[it.sample_reference];
      [
        run_id: it.scpca_run_id,
        library_id: it.scpca_library_id,
        sample_id: it.scpca_sample_id.split(";").sort().join(","),
        project_id: Utils.parseNA(it.scpca_project_id)?: "no_project",
        submitter: Utils.parseNA(it.submitter),
        technology: it.technology,
        assay_ontology_term_id: Utils.parseNA(it.assay_ontology_term_id),
        seq_unit: it.seq_unit,
        submitter_cell_types_file: Utils.parseNA(it.submitter_cell_types_file),
        openscpca_cell_types_file: Utils.parseNA(it.openscpca_cell_types_file),
        feature_barcode_file: Utils.parseNA(it.feature_barcode_file),
        feature_barcode_geom: Utils.parseNA(it.feature_barcode_geom),
        files_directory: Utils.parseNA(it.files_directory),
        slide_serial_number: Utils.parseNA(it.slide_serial_number),
        slide_section: Utils.parseNA(it.slide_section),
        ref_assembly: it.sample_reference,
        ref_fasta: "${params.ref_rootdir}/${sample_refs.ref_fasta}",
        ref_fasta_index: "${params.ref_rootdir}/${sample_refs.ref_fasta_index}",
        ref_gtf: "${params.ref_rootdir}/${sample_refs.ref_gtf}",
        mito_file: "${params.ref_rootdir}/${sample_refs.mito_file}",
        // account for the refs sometimes being null
        salmon_splici_index: sample_refs.splici_index ? "${params.ref_rootdir}/${sample_refs.splici_index}" : '',
        t2g_3col_path: sample_refs.t2g_3col_path ? "${params.ref_rootdir}/${sample_refs.t2g_3col_path}" : '',
        salmon_bulk_index: sample_refs.salmon_bulk_index ? "${params.ref_rootdir}/${sample_refs.salmon_bulk_index}" : '',
        t2g_bulk_path: sample_refs.t2g_bulk_path ? "${params.ref_rootdir}/${sample_refs.t2g_bulk_path}" : '',
        cellranger_index: sample_refs.cellranger_index ? "${params.ref_rootdir}/${sample_refs.cellranger_index}" : '',
        star_index: sample_refs.star_index ? "${params.ref_rootdir}/${sample_refs.star_index}" : '',
        infercnv_gene_order: sample_refs.infercnv_gene_order ? "${params.ref_rootdir}/${sample_refs.infercnv_gene_order}" : '',
        scpca_version: workflow.revision ?: workflow.manifest.version,
        nextflow_version: nextflow.version.toString()
      ]
    }

  // branch runs by technology
  runs_ch = unfiltered_runs_ch
    .branch{
      bulk: it.technology in bulk_techs
      feature: (it.technology in citeseq_techs) || (it.technology in cellhash_techs)
      rna: it.technology in rna_techs
      spatial: it.technology in spatial_techs
      flex: it.technology in flex_techs
    }

  // generate lists of library ids for feature libraries & RNA-only
  feature_libs = runs_ch.feature
    .collect{it.library_id}
  rna_only_libs = runs_ch.rna
    .filter{!(it.library_id in feature_libs.getVal())}
    .collect{it.library_id}
  multiplex_libs = runs_ch.rna
    .filter{it.sample_id.contains(",")}
    .collect{it.library_id}

  // get list of samples with bulk RNA-seq
  bulk_samples = runs_ch.bulk
    .collect{it.sample_id}

  // get genetic multiplex libs with all bulk samples present
  genetic_multiplex_libs = runs_ch.rna
    .filter{!params.skip_genetic_demux} // empty channel if skipping genetic demux
    .filter{it.sample_id.contains(",")}
    .filter{it.sample_id.tokenize(",").every{it in bulk_samples.getVal()}}
    .collect{it.library_id}

  // **** Process Bulk RNA-seq data ***
  bulk_quant_rna(runs_ch.bulk)

  // **** Process Spatial Transcriptomics data ****
  spaceranger_quant(runs_ch.spatial)

  // **** Process 10x flex RNA-seq data ***
  flex_quant(runs_ch.flex, flex_probesets, file(params.cellhash_pool_file))
  flex_sce_ch = generate_sce_cellranger(flex_quant.out, file(params.sample_metafile))
    .branch{
      continue_processing: it[2].size() > 0 || it[2].name.startsWith("STUBL")
      skip_processing: true
      }

  // send library ids in flex_sce_ch.skip_processing to log
  flex_sce_ch.skip_processing
    .subscribe{
      log.error("There are no cells found in the filtered object for ${it[0].library_id}.")
    }

  // **** Process 10x tag-based RNA-seq data ****
  map_quant_rna(runs_ch.rna, cell_barcodes)

  // get RNA-only libraries
  rna_quant_ch = map_quant_rna.out
    .filter{it[0]["library_id"] in rna_only_libs.getVal()}
  // make rds for rna only
  rna_sce_ch = generate_sce(rna_quant_ch, file(params.sample_metafile))
    // only continue processing any samples with > 0 cells left after filtering
    .branch{
      continue_processing: it[2].size() > 0 || it[2].name.startsWith("STUBL")
      skip_processing: true
      }

  // send library ids in rna_sce_ch.skip_processing to log
  rna_sce_ch.skip_processing
    .subscribe{
      log.error("There are no cells found in the filtered object for ${it[0].library_id}.")
    }

  // **** Process feature data ****
  map_quant_feature(runs_ch.feature, cell_barcodes)

  // join feature & RNA quants for feature reads
  feature_rna_quant_ch = map_quant_feature.out
    .map{[it[0]["library_id"]] + it } // add library_id from metadata as first element
    // join rna quant to feature quant by library_id; expect mismatches for rna-only, so don't fail
    .join(map_quant_rna.out.map{[it[0]["library_id"]] + it },
          by: 0, failOnDuplicate: true, failOnMismatch: false)
    .map{it.drop(1)} // remove library_id index

  // make rds for RNA with feature quants
  all_feature_ch = generate_sce_with_feature(feature_rna_quant_ch, file(params.sample_metafile))
    .branch{
      continue_processing: it[2].size() > 0 || it[2].name.startsWith("STUB")
      skip_processing: true
    }

  // send library ids in all_feature_ch.skip_processing to log
  all_feature_ch.skip_processing
    .subscribe{
      log.error("There are no cells found in the filtered object for ${it[0].library_id}.")
    }

  // pull out cell hash libraries for demuxing
  feature_sce_ch = all_feature_ch.continue_processing
    .branch{ // branch cellhash libs
      cellhash: it[0]["feature_meta"]["technology"] in cellhash_techs
      single: true
    }

  // apply cellhash demultiplexing
  cellhash_demux_ch = cellhash_demux_sce(feature_sce_ch.cellhash, file(params.cellhash_pool_file ?: empty_file))
  combined_feature_sce_ch = cellhash_demux_ch.mix(feature_sce_ch.single)

  // join SCE outputs and branch by genetic multiplexing
  sce_ch = rna_sce_ch.continue_processing.mix(combined_feature_sce_ch)
    .branch{
      genetic_multiplex: it[0]["library_id"] in genetic_multiplex_libs.getVal()
      no_genetic: true
    }

  // **** Perform Genetic Demultiplexing ****
  genetic_multiplex_run_ch = runs_ch.rna
    .filter{it.library_id in genetic_multiplex_libs.getVal()}
  genetic_demux_vireo(genetic_multiplex_run_ch, unfiltered_runs_ch, cell_barcodes, bulk_techs)


  // join demux result with SCE output (fail if there are any missing or extra libraries)
  // output structure: [meta_demux, vireo_dir, meta_sce, sce_rds]
  demux_results_ch = genetic_demux_vireo.out
    .map{[it[0]["library_id"]] + it }
    .join(sce_ch.genetic_multiplex.map{[it[0]["library_id"]] + it },
          by: 0, failOnDuplicate: true, failOnMismatch: true)
    .map{it.drop(1)}
  // add genetic demux results to sce objects
  genetic_demux_sce(demux_results_ch)

  // **** Post processing and generate QC reports ****
  // combine all SCE outputs
  // Make channel for all library sce files
  all_sce_ch = sce_ch.no_genetic.mix(flex_sce_ch.continue_processing, genetic_demux_sce.out)
    // add  unique ID  to metadata that will be used to label output folders and join skipped libraries throughout workflow
    .map{ meta_in, unfiltered_sce, filtered_sce ->
      def meta = meta_in.clone(); // clone meta before adding in unique id
      // we can't use library ID for flex multiplexed so will use sample and library ID for those samples only
      meta.unique_id = (meta.technology in ["10Xflex_v1.1_multi"]) ? "${meta.library_id}-${meta.sample_id}" : "${meta.library_id}";
      // return updated meta and sce files
      return [
        meta,
        unfiltered_sce,
        filtered_sce
       ]
    }
  post_process_sce(all_sce_ch)


  post_process_ch = post_process_sce.out
    // only continue processing any samples with > 0 cells left after processing
    .branch{
      continue_processing: it[3].size() > 0 || it[3].name.startsWith("STUB")
      skip_processing: true
      }

  // send library ids in post_process_ch.skip_processing to log
  post_process_ch.skip_processing
    .subscribe{
      log.error("There are no cells found in the processed object for ${it[0].library_id}.")
    }

  // Cluster SCE
  cluster_sce(post_process_ch.continue_processing)

 // Perform celltyping and call CNVs, if specified
  if (perform_celltyping) { // use perform_celltyping, not params.perform_celltyping
    annotated_celltype_ch = annotate_celltypes(cluster_sce.out)
    if (params.perform_cnv_inference) {
     cnv_celltype_ch = call_cnvs(annotated_celltype_ch)
    }
  } else {
    annotated_celltype_ch = cluster_sce.out
  }


  // If CNV inference not requested, make sure there's an empty heatmap placeholder
  if (!params.perform_cnv_inference) {
    cnv_celltype_ch = annotated_celltype_ch
      .map{ meta, unfiltered, filtered, processed -> tuple(
        meta,
        unfiltered,
        filtered,
        processed,
        empty_file)
    }
  }

  cnv_celltype_ch.view()

  // first mix any skipped libraries from both rna and feature libs
  no_filtered_ch = rna_sce_ch.skip_processing.mix(all_feature_ch.skip_processing, flex_sce_ch.skip_processing)
    // add a fake processed file
    .map{meta, unfiltered, filtered -> tuple(
      meta,
      unfiltered,
      filtered,
      empty_file,
      [] // infercnv heatmap empty file placeholder, but don't repeat empty_file
    )}

  // combine back with libraries that skipped filtering and post processing
  sce_output_ch = post_process_ch.skip_processing
    .map{ meta, unfiltered, filtered, processed -> tuple(
        meta,
        unfiltered_sce,
        filtered_sce,
        processed_sce,
        empty_file)
    }
    .mix(cnv_celltype_ch, no_filtered_ch)

  def report_template_tuple = tuple(
    file(params.report_template_dir, type: 'dir', checkIfExists: true),
    params.report_template_file,
    params.celltype_report_template_file
  )

  // generate QC reports & metrics, then publish sce
  qc_publish_sce(
    sce_output_ch,
    report_template_tuple,
    // need to pass this in as a value since the param may have been overridden
    perform_celltyping,
    // paths to files needed to make consensus cell type validation dot plots
    file(params.validation_groups_file),
    file(params.validation_markers_file),
    file(params.validation_palette_file)
  )

  // convert SCE object to anndata
  anndata_ch = qc_publish_sce.out.data
    // skip multiplexed libraries
    .filter{!(it[0]["library_id"] in multiplex_libs.getVal())}
  sce_to_anndata(anndata_ch)

}
