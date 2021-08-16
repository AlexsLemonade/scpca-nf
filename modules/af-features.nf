


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
  tag "${run_id}-features"
  input:
    tuple val(sample_id), val(run_id), val(tech), 
          path(read1), path(read2), 
          val(feature_geom), path(feature_index)
  output:
    tuple val(sample_id), val(run_id), 
          path(run_dir), path(feature_index)
  script:
    // label the run directory by id
    run_dir = "${run_id}-features"
    // Define umi geometries
    umi_geoms = ['CITEseq_10Xv2': '1[17-26]',
                 'CITEseq_10Xv3': '1[17-28]',
                 'CITEseq_10Xv3.1': '1[17-28]']
    """
    mkdir -p ${run_dir}
    salmon alevin \
      -l ISR \
      -1 ${read1} \
      -2 ${read2} \
      -i ${feature_index} \
      --read-geometry ${feature_geom} \
      --bc-geometry 1[1-16] \
      --umi-geometry ${umi_geoms[tech]} \
      --rad \
      -o ${run_dir} \
      -p ${task.cpus} 
    """
}

// quantify features from rad input
process fry_quant_feature{
  container params.ALEVINFRY_CONTAINER
  label 'cpus_8'
  publishDir "${params.outdir}/features"

  input:
    tuple val(sample_id), val(run_id),
          path(run_dir), path(feature_index)
    path barcode_file
  output:
    tuple val(sample_id), val(run_id),
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
      --tg-map ${feature_index}/t2g.tsv \
      --resolution ${params.resolution} \
      -o ${run_dir} \
      --use-mtx \
      -t ${task.cpus} \

    # remove large files
    rm ${run_dir}/*.rad ${run_dir}/*.bin 
    """
}


workflow map_quant_feature{
  take: feature_channel
  // a channel with a groovy map for each feature barcode library to process
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
      .map{file("${params.barcode_dir}/${params.cell_barcodes[it.technology]}")}

    // run Alevin on feature reads
    alevin_feature(feature_reads_ch)
    // quantify feature reads 
    fry_quant_feature(alevin_feature.out, feature_cellbarcode_ch)
  
  emit: fry_quant_feature.out
  // a tuple of sample_id, run_id, and the alevin-fry output directory
}
