
//index a feature barcode file
process index_feature{
  container params.SALMON_CONTAINER
  
  input:
    tuple val(id), path(feature_file)
  output:
    tuple val(id), path("feature_index")
  script:
    """
    salmon index \
      -t ${feature_file} \
      -i feature_index \
      --features \
      -k 7 

    awk '{print \$1"\\t"\$1;}' ${feature_file} > feature_index/t2g.tsv
    """
}

// generates RAD file for alevin feature matrix using alevin
process alevin_feature{
  container params.SALMON_CONTAINER
  label 'cpus_8'
  label 'mem_8'
  tag "${meta.run_id}-features"
  publishDir "${params.checkpoints_dir}/rad/${meta.library_id}", enabled: params.publish_fry_outs
  input:
    tuple val(meta), 
          path(read1), path(read2), 
          path(feature_index)
  output:
    tuple val(meta),
          path(run_dir)
  script:
    // label the run directory by id
    run_dir = "${meta.run_id}-features"
    // Define umi geometry by 10x version
    umi_geom_map = ['10Xv2': '1[17-26]',
                    '10Xv3': '1[17-28]',
                    '10Xv3.1': '1[17-28]']
    tech_version = meta.technology.split('_').last()
    umi_geom = umi_geom_map[tech_version]
    """
    mkdir -p ${run_dir}
    salmon alevin \
      -l ISR \
      -1 ${read1} \
      -2 ${read2} \
      -i ${feature_index} \
      --read-geometry ${meta.feature_barcode_geom} \
      --bc-geometry 1[1-16] \
      --umi-geometry ${umi_geom} \
      --rad \
      -o ${run_dir} \
      -p ${task.cpus}

    cp ${feature_index}/t2g.tsv ${run_dir}/t2g.tsv
    """
}

// quantify features from rad input
process fry_quant_feature{
  container params.ALEVINFRY_CONTAINER
  label 'cpus_8'
  tag "${meta.run_id}-features"
  publishDir "${params.checkpoints_dir}/af/${meta.library_id}", enabled: params.publish_fry_outs
  input:
    tuple val(meta),
          path(run_dir)
    path barcode_file
  output:
    tuple val(meta),
          path(run_dir)
  
  script: 
    """
    alevin-fry generate-permit-list \
      -i ${run_dir} \
      --expected-ori fw \
      -o ${run_dir} \
      --unfiltered-pl ${barcode_file}

    alevin-fry collate \
      --input-dir ${run_dir} \
      --rad-dir ${run_dir} \
      -t ${task.cpus}
    
    alevin-fry quant \
      --input-dir ${run_dir} \
      --tg-map ${run_dir}/t2g.tsv \
      --resolution ${params.af_resolution} \
      -o ${run_dir} \
      --use-mtx \
      -t ${task.cpus} \

    # remove large files
    rm ${run_dir}/*.rad ${run_dir}/*.bin 
    """
}


workflow map_quant_feature{
  take: feature_channel
  // a channel with a groovy map of metadata for each feature barcode library to process
  main:
    //get and map the feature barcode files
    feature_barcodes_ch = feature_channel
      .map{meta -> tuple(meta.feature_barcode_file,
                         file("${meta.feature_barcode_file}"))}
      .unique()
    index_feature(feature_barcodes_ch)

    // create tuple of [metadata, [Read1 files], [Read2 files]]
    // We start by including the feature_barcode file so we can join to the indices, but that will be removed
    feature_reads_ch = feature_channel
      .map{meta -> tuple(meta.feature_barcode_file,
                         meta,
                         file("${meta.files_directory}/*_R1_*.fastq.gz"),
                         file("${meta.files_directory}/*_R2_*.fastq.gz")
                        )}
      .combine(index_feature.out, by: 0) // combine by the feature_barcode_file
      .map{ it.subList(1, it.size())} // remove the first element (feature_barcode_file)
    
    cellbarcode_ch = feature_channel
      .map{file("${params.barcode_dir}/${params.cell_barcodes[it.technology]}")}

    // run Alevin on feature reads
    alevin_feature(feature_reads_ch)
    // quantify feature reads 
    fry_quant_feature(alevin_feature.out, cellbarcode_ch)
  
  emit: fry_quant_feature.out
  // a tuple of metadata map and the alevin-fry output directory
}
