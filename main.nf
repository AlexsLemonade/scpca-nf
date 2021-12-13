#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// run_ids are comma separated list to be parsed into 
// a list of run ids, library ids, and or sample_ids
// or "All" to process all samples in the metadata file
params.run_ids = "SCPCR000001,SCPCS000101,SCPCR000050,SCPCR000084"
// to run all samples in a project use that project's submitter name
params.project = ""

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
single_cell_techs= cell_barcodes.keySet()
bulk_techs = ['single_end', 'paired_end']
spatial_techs = ["spatial", "visium_v1", "visium_v2"]
all_techs = single_cell_techs + bulk_techs + spatial_techs
rna_techs = single_cell_techs.findAll{it.startsWith('10Xv')}
feature_techs = single_cell_techs.findAll{it.startsWith('CITEseq') || it.startsWith('cellhash')}


// include processes from modules
include { map_quant_rna } from './modules/af-rna.nf' addParams(cell_barcodes: cell_barcodes)
include { map_quant_feature } from './modules/af-features.nf' addParams(cell_barcodes: cell_barcodes)
include { bulk_quant_rna } from './modules/bulk-salmon.nf'
include { spaceranger_quant } from './modules/spaceranger.nf'
include { generate_sce; generate_merged_sce } from './modules/generate-rds.nf'
include { sce_qc_report } from './modules/qc-report.nf'

workflow {
  // select runs to use
  if (params.project){
    // projects will use all runs in the project & supersede run_ids
    run_ids = []
  }else{
    run_ids = params.run_ids?.tokenize(',') ?: []
  }
  run_all = run_ids[0] == "All"

  runs_ch = Channel.fromPath(params.run_metafile)
    .splitCsv(header: true, sep: '\t')
    // convert row data to a metadata map, keeping only columns we will need (& some renaming)
    .map{[
      run_id: it.scpca_run_id,
      library_id: it.scpca_library_id,
      sample_id: it.scpca_sample_id,
      project_id: it.scpca_project_id?: "no_project",
      submitter: it.submitter,
      technology: it.technology,
      seq_unit: it.seq_unit,
      feature_barcode_file: it.feature_barcode_file,
      feature_barcode_geom: it.feature_barcode_geom,
      s3_prefix: it.s3_prefix,
      slide_serial_number: it.slide_serial_number,
      slide_section: it.slide_section,
      files: it.files
    ]}
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
  
  // generate lists of library ids for feature libraries & RNA-only
  feature_libs = runs_ch.filter{it.technology in feature_techs}
    .collect{it.library_id}
  rna_only_libs = runs_ch.filter{!(it.library_id in feature_libs.getVal())}
    .collect{it.library_id}

  // **** Process Bulk RNA-seq data *** 
  bulk_ch = runs_ch.filter{it.technology in bulk_techs}
  bulk_quant_rna(bulk_ch)

  // **** Process RNA-seq data ****
  rna_ch = runs_ch.filter{it.technology in rna_techs}
  map_quant_rna(rna_ch)
  
  
  // get RNA-only libraries
  rna_quant_ch = map_quant_rna.out
    .filter{it[0]["library_id"] in rna_only_libs.getVal()}
  // make rds for rna only
  rna_sce_ch = generate_sce(rna_quant_ch)


  // **** Process feature data ****
  feature_ch = runs_ch.filter{it.technology in feature_techs} 
  map_quant_feature(feature_ch)
  
  // combine feature & RNA quants for feature reads
  feature_rna_quant_ch = map_quant_feature.out
    .map{[it[0]["library_id"]] + it } // add library_id from metadata as first element
    .combine(map_quant_rna.out.map{[it[0]["library_id"]] + it }, by: 0) // combine by library_id 
    .map{it.subList(1, it.size())} // remove library_id index
  // make rds for merged RNA and feature quants
  merged_sce_ch = generate_merged_sce(feature_rna_quant_ch)

   // **** Process Spatial Transcriptomics data ****
  spatial_ch = runs_ch.filter{it.technology in spatial_techs}
  spaceranger_quant(spatial_ch)

  // **** Generate QC reports ****
  // Make channel for all library sce files & run QC report
  library_sce_ch = rna_sce_ch.mix(merged_sce_ch)
  sce_qc_report(library_sce_ch)
}
