// containers
ALEVINFRY_CONTAINER = 'quay.io/biocontainers/alevin-fry:0.4.1--h7d875b9_0'

// quantify rna from rad input
process fry_quant_rna{
  container ALEVINFRY_CONTAINER
  label 'cpus_8'
  publishDir "${params.outdir}/rna"

  input:
    tuple val(sample_id), val(run_id), path(run_dir)
    path barcode_file
    path tx2gene_3col
  output:
    tuple val(sample_id), val(run_id), path(run_dir)
  
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
      --tg-map ${tx2gene_3col} \
      --resolution ${params.resolution} \
      -o ${run_dir} \
      --use-mtx \
      -t ${task.cpus} \

    # remove large files
    rm ${run_dir}/*.rad ${run_dir}/*.bin 
    """
}




// quantify features from rad input
process fry_quant_feature{
  container ALEVINFRY_CONTAINER
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
