#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// 10X barcode files
cell_barcodes = [
  '10Xv2': '737K-august-2016.txt',
  '10Xv2_5prime': '737K-august-2016.txt',
  '10Xv3': '3M-february-2018.txt',
  '10Xv3.1': '3M-february-2018.txt',
  'CITEseq_10Xv2': '737K-august-2016.txt',
  'CITEseq_10Xv3': '3M-february-2018.txt',
  'CITEseq_10Xv3.1': '3M-february-2018.txt',
  'cellhash_10Xv2': '737K-august-2016.txt',
  'cellhash_10Xv3': '3M-february-2018.txt',
  'cellhash_10Xv3.1': '3M-february-2018.txt'
  ]

// supported technologies
single_cell_techs = cell_barcodes.keySet()
bulk_techs = ['single_end', 'paired_end']
spatial_techs = ["spatial", "visium_v1", "visium_v2"]
all_techs = single_cell_techs + bulk_techs + spatial_techs
rna_techs = single_cell_techs.findAll{it.startsWith('10Xv')}
feature_techs = single_cell_techs.findAll{it.startsWith('CITEseq') || it.startsWith('cellhash')}


// include processes from modules
include { map_quant_rna } from './modules/af-rna.nf' addParams(cell_barcodes: cell_barcodes)
include { map_quant_feature } from './modules/af-features.nf' addParams(cell_barcodes: cell_barcodes)
include { bulk_quant_rna } from './modules/bulk-salmon.nf'
include { genetic_demux } from './modules/genetic-demux.nf' addParams(cell_barcodes: cell_barcodes, bulk_techs: bulk_techs)
include { spaceranger_quant } from './modules/spaceranger.nf'
include { generate_sce; generate_merged_sce; multiplex_demux_sce } from './modules/sce-processing.nf'
include { sce_qc_report } from './modules/qc-report.nf'

// parameter checks
param_error = false

if (!file(params.run_metafile).exists()) {
  log.error("The 'run_metafile' file '${params.run_metafile}' can not be found.")
  param_error = true
}

resolution_strategies = ['cr-like', 'full', 'cr-like-em', 'parsimony', 'trivial']
if (!resolution_strategies.contains(params.af_resolution)) {
  log.error("'af_resolution' must be one of the following: ${resolution_strategies}")
  param_error = true
}

if(param_error){
  System.exit(1)
}

// Main workflow
workflow {
  // select runs to use
  if (params.project){
    // projects will use all runs in the project & supersede run_ids
    run_ids = []
  }else{
    run_ids = params.run_ids?.tokenize(',') ?: []
  }
  run_all = run_ids[0] == "All"
  if (run_all){
    log.info("Executing workflow for all runs in the run metafile.")
  }

  unfiltered_runs_ch = Channel.fromPath(params.run_metafile)
    .splitCsv(header: true, sep: '\t')
    // convert row data to a metadata map, keeping only columns we will need (& some renaming)
    .map{[
      run_id: it.scpca_run_id,
      library_id: it.scpca_library_id,
      sample_id: it.scpca_sample_id.split(";").sort().join(","),
      project_id: it.scpca_project_id?: "no_project",
      submitter: it.submitter,
      technology: it.technology,
      seq_unit: it.seq_unit,
      feature_barcode_file: it.feature_barcode_file,
      feature_barcode_geom: it.feature_barcode_geom,
      files_directory: it.files_directory,
      slide_serial_number: it.slide_serial_number,
      slide_section: it.slide_section,
      files: it.files
    ]}

 runs_ch = unfiltered_runs_ch 
    // only technologies we know how to process
    .filter{it.technology in all_techs} 
    // use only the rows in the run_id list (run, library, or sample can match)
    // or run by project or submitter if the project parameter is set
    .filter{run_all 
             || (it.run_id in run_ids) 
             || (it.library_id in run_ids)
             || (it.sample_id in run_ids)
             || (it.submitter == params.project)
             || (it.project_id == params.project)
            }
     .branch{
       bulk: it.technology in bulk_techs
       feature: it.technology in feature_techs
       rna: it.technology in rna_techs 
       spatial: it.technology in spatial_techs
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

  // **** Process Bulk RNA-seq data *** 
  bulk_quant_rna(runs_ch.bulk)

  // **** Process RNA-seq data ****
  map_quant_rna(runs_ch.rna)
  
  // get RNA-only libraries
  rna_quant_ch = map_quant_rna.out
    .filter{it[0]["library_id"] in rna_only_libs.getVal()}
  // make rds for rna only
  rna_sce_ch = generate_sce(rna_quant_ch)


  // **** Process feature data ****
  map_quant_feature(runs_ch.feature)
  
  // combine feature & RNA quants for feature reads
  feature_rna_quant_ch = map_quant_feature.out
    .map{[it[0]["library_id"]] + it } // add library_id from metadata as first element
    .combine(map_quant_rna.out.map{[it[0]["library_id"]] + it }, by: 0) // combine by library_id 
    .map{it.subList(1, it.size())} // remove library_id index
  // make rds for merged RNA and feature quants
  merged_sce_ch = generate_merged_sce(feature_rna_quant_ch)

  // join SCE outputs and branch by multiplexing
  sce_ch = rna_sce_ch.mix(merged_sce_ch)
    .branch{
      multiplex: it[0]["library_id"] in multiplex_libs.getVal()
      single: true
    }

  // **** Process multiplexed samples ****
  multiplex_run_ch = runs_ch.rna
    .filter{it.library_id in multiplex_libs.getVal()} 
  genetic_demux(multiplex_run_ch, unfiltered_runs_ch, sce_ch.multiplex)
  // combine demux result with SCE output
  // output structure: [meta_demux, vireo_dir, meta_sce, sce_rds]
  sce_demux_ch = genetic_demux.out
    .map{[it[0]["library_id"]] + it }
    .combine(sce_ch.multiplex.map{[it[0]["library_id"]] + it }, by: 0)
    .map{it.subList(1, it.size())}
  multiplex_demux_sce(sce_demux_ch, file(params.cellhash_pool_file))

  // **** Generate QC reports ****
  // combine all SCE outputs
  // Make channel for all library sce files & run QC report
  all_sce_ch = sce_ch.single.mix(multiplex_demux_sce.out)
  sce_qc_report(all_sce_ch)


   // **** Process Spatial Transcriptomics data ****
  spaceranger_quant(runs_ch.spatial)
}
