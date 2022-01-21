

// generates RAD file using alevin
process alevin_rad{
  container params.SALMON_CONTAINER
  label 'cpus_12'
  label 'disk_dynamic'
  tag "${meta.run_id}-rna"
  publishDir "${params.outdir}/internal/rad/${meta.library_id}"
  input:
    tuple val(meta), 
          path(read1), path(read2)
    path index
  output:
    tuple val(meta), path(run_dir) 
  script:
    // label the run-dir
    run_dir = "${meta.run_id}-rna"
    // choose flag by technology
    tech_flag = ['10Xv2': '--chromium',
                 '10Xv2_5prime': '--chromium',
                 '10Xv3': '--chromiumV3',
                 '10Xv3.1': '--chromiumV3']
    // run alevin like normal with the --rad flag 
    // creates output directory with RAD file needed for alevin-fry
    """
    mkdir -p ${run_dir}
    salmon alevin \
      -l ISR \
      ${tech_flag[meta.technology]} \
      -1 ${read1} \
      -2 ${read2} \
      -i ${index} \
      -o ${run_dir} \
      -p ${task.cpus} \
      --dumpFeatures \
      --rad
    """
}

// quantify rna from RAD input
process fry_quant_rna{
  container params.ALEVINFRY_CONTAINER
  label 'cpus_8'
  tag "${meta.run_id}-rna"
  publishDir "${params.outdir}/internal/af/${meta.library_id}"

  input:
    tuple val(meta), path(run_dir), path(barcode_file)
    path tx2gene_3col
  output:
    tuple val(meta), path(run_dir)
  
  script: 
    """
    alevin-fry generate-permit-list \
      -i ${run_dir} \
      --expected-ori ${meta.technology == '10Xv2_5prime' ? 'rc' : 'fw'} \
      -o ${run_dir} \
      --unfiltered-pl ${barcode_file}

    alevin-fry collate \
      --input-dir ${run_dir} \
      --rad-dir ${run_dir} \
      -t ${task.cpus}
    
    alevin-fry quant \
      --input-dir ${run_dir} \
      --tg-map ${tx2gene_3col} \
      --resolution ${params.af_resolution} \
      -o ${run_dir} \
      --use-mtx \
      -t ${task.cpus} \

    # remove large files
    rm ${run_dir}/*.rad ${run_dir}/*.bin 
    """
}


workflow map_quant_rna {
  take: rna_channel
  // a channel with a map of metadata for each rna library to process
  main:
    // create tuple of (metadata map, [Read1 files], [Read2 files])
    // for rnaseq runs
    rna_channel = rna_channel
      .map{it.rad_publish_dir = "${params.outdir}/internal/rad/${it.library_id}";
           it.rad_dir = "${it.rad_publish_dir}/${it.run_id}-rna"; 
           it.barcode_file = "${params.barcode_dir}/${params.cell_barcodes[it.technology]}";
           it}

    rna_reads_ch = rna_channel
      .filter{!(params.rad_skip && file(it.rad_dir).exists())}      
      .map{meta -> tuple(meta,
                         file("s3://${meta.s3_prefix}/*_R1_*.fastq.gz"),
                         file("s3://${meta.s3_prefix}/*_R2_*.fastq.gz")
                        )}

    rna_rad_ch = rna_channel
      .filter{params.rad_skip && file(it.rad_dir).exists()}
      .map{meta -> tuple(meta, 
                         file(meta.rad_dir)
                         )}

    // run Alevin for mapping
    alevin_rad(rna_reads_ch, params.splici_index)
    // quantify with alevin-fry
    all_rad_ch = alevin_rad.out.mix(rna_rad_ch)
      .map{it + [file(it[0].barcode_file)]}

    fry_quant_rna(all_rad_ch, params.t2g_3col_path)
  
  emit: fry_quant_rna.out
  // a tuple of meta and the alevin-fry output directory
}
