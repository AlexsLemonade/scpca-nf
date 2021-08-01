

// containers
SALMON_CONTAINER = 'quay.io/biocontainers/salmon:1.5.2--h84f40af_0'

// generates RAD file using alevin
process alevin_rad{
  container SALMON_CONTAINER
  label 'cpus_12'
  tag "${run_id}-rna"
  input:
    tuple val(sample_id), val(run_id), val(tech), 
          path(read1), path(read2)
    path index
  output:
    tuple val(sample_id), val(run_id), path(run_dir)
  script:
    // label the run-dir
    run_dir = "${run_id}-rna"
    // choose flag by technology
    tech_flag = ['10Xv2': '--chromium',
                 '10Xv3': '--chromiumV3',
                 '10Xv3.1': '--chromiumV3']
    // run alevin like normal with the --rad flag 
    // creates output directory with RAD file needed for alevin-fry
    """
    mkdir -p ${run_dir}
    salmon alevin \
      -l ISR \
      ${tech_flag[tech]} \
      -1 ${read1} \
      -2 ${read2} \
      -i ${index} \
      -o ${run_dir} \
      -p ${task.cpus} \
      --dumpFeatures \
      --rad
    """
}

// generates RAD file for alevin feature matrix using alevin
process alevin_feature{
  container SALMON_CONTAINER
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

process index_feature{
  container SALMON_CONTAINER
  
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
