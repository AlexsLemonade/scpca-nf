#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// run parameters
params.barcode_dir = 's3://nextflow-ccdl-data/reference/10X/barcodes' 
params.run_metafile = 's3://ccdl-scpca-data/sample_info/scpca-library-metadata.tsv'
params.outdir = "s3://nextflow-ccdl-results/scpca/processed"

// file paths
params.index_path = 's3://nextflow-ccdl-data/reference/homo_sapiens/ensembl-103/salmon_index/spliced_intron_txome_k31'
params.t2g_2col_path = 's3://nextflow-ccdl-data/reference/homo_sapiens/ensembl-103/annotation/Homo_sapiens.GRCh38.103.spliced_intron.tx2gene.tsv'
params.t2g_3col_path = 's3://nextflow-ccdl-data/reference/homo_sapiens/ensembl-103/annotation/Homo_sapiens.GRCh38.103.spliced_intron.tx2gene_3col.tsv'

// run_ids are comma separated list to be parsed into a list of run ids,
// or "All" to process all samples in the metadata file
params.run_ids = "SCPCR000001,SCPCR000002"

// 10X barcode files
barcodes = ['10Xv2': '737K-august-2016.txt',
            '10Xv3': '3M-february-2018.txt',
            '10Xv3.1': '3M-february-2018.txt']

// supported single cell technologies
tech_list = barcodes.keySet()

// generates RAD file using alevin
process alevin{
  container 'quay.io/biocontainers/salmon:1.5.1--h84f40af_0'
  label 'cpus_8'
  tag "${id}"
  publishDir "${params.outdir}"
  input:
    tuple val(id), val(tech), path(read1), path(read2)
    path index
    path tx2gene_2col
  output:
    path run_dir
  script:
    // label the run-dir
    run_dir = "${id}"
    // choose flag by technology
    tech_flag = ['10Xv2': '--chromium',
                 '10Xv3': '--chromiumV3',
                 '10Xv3.1': '--chromiumV3']
    // run alevin like normal with the --justAlign flag 
    // creates output directory with RAD file needed for alevin-fry
    """
    mkdir -p ${run_dir}
    salmon alevin \
      -l ISR \
      ${tech_flag[tech]} \
      -1 ${read1} \
      -2 ${read2} \
      -i ${index} \
      --tgMap ${tx2gene_2col} \
      -o ${run_dir} \
      -p ${task.cpus} \
      --dumpFeatures \
      --justAlign
    """
}

//generate permit list from RAD input 
process generate_permit{
  container 'quay.io/biocontainers/alevin-fry:0.4.0--h7d875b9_0'
  input:
    path run_dir
    path barcode_file
  output:
    path run_dir
  script:
    """
    alevin-fry generate-permit-list \
      -i ${run_dir} \
      --expected-ori fw \
      -o ${run_dir} \
      --unfiltered-pl ${barcode_file}
    """
}

// given permit list and barcode mapping, collate RAD file 
process collate_fry{
  container 'quay.io/biocontainers/alevin-fry:0.4.0--h7d875b9_0'
  label 'cpus_8'
  input: 
    path run_dir
  output: 
    path run_dir
  script:
    """
    alevin-fry collate \
      --input-dir ${run_dir} \
      --rad-dir ${run_dir} \
      -t ${task.cpus}
    rm ${run_dir}/map.rad
    """
}

// then quantify collated RAD file
process quant_fry{
  container 'quay.io/biocontainers/alevin-fry:0.4.0--h7d875b9_0'
  label 'cpus_8'
  publishDir "${params.outdir}"
  input: 
    path run_dir
    path tx2gene_3col
  output: 
    path run_dir
  script:
    """
    alevin-fry quant \
     --input-dir ${run_dir} \
     --tg-map ${tx2gene_3col} \
     --output-dir ${run_dir} \
     --resolution cr-like \
     -t ${task.cpus} \
     --use-mtx
    rm ${run_dir}/*.rad
    rm ${run_dir}/*.bin
    """
}

// insert process for generating unfiltered and filtered output in desired format 
// insert process for generating QC report 

workflow{
  run_ids = params.run_ids?.tokenize(',') ?: []
  run_all = run_ids[0] == "All"
  samples_ch = Channel.fromPath(params.run_metafile)
    .splitCsv(header: true, sep: '\t')
    .filter{it.technology in tech_list} 
    // use only the rows in the sample list
    .filter{run_all || (it.scpca_run_id in run_ids)}
  // create tuple of [sample_id, technology, [Read1 files], [Read2 files]]
  reads_ch = samples_ch
    .map{row -> tuple(row.scpca_run_id,
                      row.technology,
                      file("s3://${row.s3_prefix}/*_R1_*.fastq.gz"),
                      file("s3://${row.s3_prefix}/*_R2_*.fastq.gz"),
                      )}

  barcodes_ch = samples_ch
    .map{row -> file("${params.barcode_dir}/${barcodes[row.technology]}")}

  // run Alevin
  alevin(reads_ch, params.index_path, params.t2g_2col_path)
  // generate permit list from alignment 
  generate_permit(alevin.out, barcodes_ch)
  // collate RAD files 
  collate_fry(generate_permit.out)
  // create gene x cell matrix
  quant_fry(collate_fry.out, params.t2g_3col_path)
}
