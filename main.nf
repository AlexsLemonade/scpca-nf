#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// run parameters
params.run_metafile = 's3://ccdl-scpca-data/sample_info/scpca-library-metadata.tsv'
params.outdir = "s3://nextflow-ccdl-results/scpca/processed"

params.resolution = 'cr-like' //default resolution is cr-like, can also use full, cr-like-em, parsimony, and trivial
params.barcode_dir = 's3://nextflow-ccdl-data/reference/10X/barcodes' 
params.index_path = 's3://nextflow-ccdl-data/reference/homo_sapiens/ensembl-103/salmon_index/spliced_intron_txome_k31'
params.t2g_3col_path = 's3://nextflow-ccdl-data/reference/homo_sapiens/ensembl-103/annotation/Homo_sapiens.GRCh38.103.spliced_intron.tx2gene_3col.tsv'

// run_ids are comma separated list to be parsed into a list of run ids,
// or "All" to process all samples in the metadata file
params.run_ids = "SCPCR000001,SCPCS000101,SCPCR000050,SCPCR000084"

// include processes from modules
include { alevin_rad; alevin_feature; index_feature } from './modules/salmon.nf'
include { fry_quant_rna; fry_quant_feature } from './modules/alevin-fry.nf'

// 10X barcode files
cell_barcodes = [
  '10Xv2': '737K-august-2016.txt',
  '10Xv3': '3M-february-2018.txt',
  '10Xv3.1': '3M-february-2018.txt',
  'CITEseq_10Xv2': '737K-august-2016.txt',
  'CITEseq_10Xv3': '3M-february-2018.txt',
  'CITEseq_10Xv3.1': '3M-february-2018.txt'
  ]

// supported single cell technologies
tech_list = cell_barcodes.keySet()
rna_techs = tech_list.findAll{it.startsWith('10Xv')}
feature_techs = tech_list.findAll{it.startsWith('CITEseq')}

workflow map_quant_rna {
  take: rna_channel
  main:
    // create tuple of [run_id, sample_id, technology, [Read1 files], [Read2 files]]
    // for rnaseq runs
    rna_reads_ch = rna_channel
      .map{row -> tuple(row.scpca_sample_id,
                        row.scpca_run_id,
                        row.technology,
                        file("s3://${row.s3_prefix}/*_R1_*.fastq.gz"),
                        file("s3://${row.s3_prefix}/*_R2_*.fastq.gz"),
                        )}

    rna_cellbarcodes_ch = rna_channel
      .map{file("${params.barcode_dir}/${cell_barcodes[it.technology]}")}

    // run Alevin
    alevin_rad(rna_reads_ch, params.index_path)
    // generate permit list from alignment 
    fry_quant_rna(alevin_rad.out, rna_cellbarcodes_ch, params.t2g_3col_path)
  
  emit: fry_quant_rna.out
}

workflow map_quant_feature {
  take: feature_channel
  main:
    //get and map the feature barcode files
    feature_barcodes_ch = feature_channel
      .map{row -> tuple(row.feature_barcode_file,
                        file("s3://${row.feature_barcode_file}"))}
      .unique()
    index_feature(feature_barcodes_ch)

    // create tuple of [run_id, sample_id, technology, [Read1 files], [Read2 files], feature_geometry, feature_index]
    // We start by including the feature_barcode file so we can join to the indices, but that will be removed
    feature_reads_ch = feature_channel
      .map{row -> tuple(row.feature_barcode_file,
                        row.scpca_sample_id,
                        row.scpca_run_id,
                        row.technology,
                        file("s3://${row.s3_prefix}/*_R1_*.fastq.gz"),
                        file("s3://${row.s3_prefix}/*_R2_*.fastq.gz"),
                        row.feature_barcode_geom
                        )}
      .combine(index_feature.out, by: 0) // combine by the feature_barcode_file
      .map{ it.subList(1, it.size())} // remove the first element
    
    feature_cellbarcode_ch = feature_channel
      .map{file("${params.barcode_dir}/${cell_barcodes[it.technology]}")}

    // run Alevin on feature reads
    alevin_feature(feature_reads_ch)
    // quantify feature reads 
    fry_quant_feature(alevin_feature.out, feature_cellbarcode_ch)
  
  emit: fry_quant_feature.out
}

workflow{
  // select runs to use
  run_ids = params.run_ids?.tokenize(',') ?: []
  run_all = run_ids[0] == "All"
  runs_ch = Channel.fromPath(params.run_metafile)
    .splitCsv(header: true, sep: '\t')
    // only technologies we know how to process
    .filter{it.technology in tech_list} 
    // use only the rows in the run_id list
    .filter{run_all || (it.scpca_run_id in run_ids)}
  
  // **** Process RNA-seq data ****
  rna_ch = runs_ch.filter{it.technology in rna_techs}
  map_quant_rna(rna_ch)

  // **** Process feature data ****
  feature_ch = runs_ch.filter{it.technology in feature_techs} 
  map_quant_feature(feature_ch)
  
  // combine feature & RNA quants for feature reads
  // just print for now
  map_quant_feature.out
    .combine(map_quant_rna.out, by: 0) // combine by the sample_id
    .view()

}
