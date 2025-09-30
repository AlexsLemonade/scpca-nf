
//index a feature barcode file
process index_feature{
  container params.SALMON_CONTAINER
  tag "${id}"

  input:
    tuple val(id), path(feature_file)
  output:
    tuple val(id), path("feature_index")
  script:
    """
    # ensure salmon only gets the first 2 columns here
    cut -f 1,2 ${feature_file} > feature_barcodes.txt

    salmon index \
      -t feature_barcodes.txt \
      -i feature_index \
      --features \
      -k 7

    awk '{print \$1"\\t"\$1;}' ${feature_file} > feature_index/t2g.tsv
    """
  stub:
    """
    mkdir -p feature_index
    """
}

// generates RAD file for alevin feature matrix using alevin
process alevin_feature {
  container params.SALMON_CONTAINER
  label 'cpus_8'
  label 'mem_8'
  tag "${meta.run_id}-features"
  publishDir "${meta.feature_rad_publish_dir}", mode: 'copy'
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
    umi_geom_map = [
      '10xv2': '1[17-26]',
      '10xv3': '1[17-28]',
      '10xv3.1': '1[17-28]',
      '10xv3_5prime': '1[17-28]',
      '10xv4': '1[17-28]'
    ]
    tech_version = meta.technology.split('_').last()
    umi_geom = umi_geom_map[tech_version]
    // get meta to write as file
    // make sure workflow version strings are correct
    meta += Utils.getVersions(workflow, nextflow)
    meta_json = Utils.makeJson(meta)
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

    echo '${meta_json}' > ${run_dir}/scpca-meta.json
    """
  stub:
    run_dir = "${meta.run_id}-features"
    meta += Utils.getVersions(workflow, nextflow)
    meta_json = Utils.makeJson(meta)
    """
    mkdir -p ${run_dir}
    echo '${meta_json}' > ${run_dir}/scpca-meta.json
    """

}

// quantify features from rad input
process fry_quant_feature {
  container params.ALEVINFRY_CONTAINER
  label 'cpus_8'
  tag "${meta.run_id}-features"
  publishDir "${params.checkpoints_dir}/alevinfry/${meta.library_id}", mode: 'copy', enabled: params.publish_fry_outs
  input:
    tuple val(meta), path(rad_dir), path(barcode_file)
  output:
    tuple val(meta), path(quant_dir)
  script:
    quant_dir = rad_dir + '_quant'
    // get meta to write as file
    meta += Utils.getVersions(workflow, nextflow)
    meta_json = Utils.makeJson(meta)
    """
    alevin-fry generate-permit-list \
      -i ${rad_dir} \
      --expected-ori fw \
      -o ${quant_dir} \
      --unfiltered-pl ${barcode_file}

    alevin-fry collate \
      --input-dir ${quant_dir} \
      --rad-dir ${rad_dir} \
      -t ${task.cpus}

    alevin-fry quant \
      --input-dir ${quant_dir} \
      --tg-map ${rad_dir}/t2g.tsv \
      --resolution ${params.af_resolution} \
      -o ${quant_dir} \
      --use-mtx \
      -t ${task.cpus} \

    # copy cmd_info and aux_info to quant directory for tracking
    cp ${rad_dir}/cmd_info.json ${quant_dir}/cmd_info.json
    cp -r ${rad_dir}/aux_info ${quant_dir}/aux_info

    echo '${meta_json}' > ${quant_dir}/scpca-meta.json
    """
  stub:
    quant_dir = rad_dir + '_quant'
    meta += Utils.getVersions(workflow, nextflow)
    meta_json = Utils.makeJson(meta)
    """
    mkdir -p '${quant_dir}'
    echo '${meta_json}' > ${quant_dir}/scpca-meta.json
    """
}


workflow map_quant_feature {
  take:
    feature_channel // a channel with a groovy map of metadata for each feature barcode library to process
    cell_barcodes // map of cell barcode files for each technology
  main:
    //get and map the feature barcode files
    feature_barcodes_ch = feature_channel
      .map{meta -> tuple(
        meta.feature_barcode_file,
        file("${meta.feature_barcode_file}", checkIfExists: true)
      )}
      .unique()
    index_feature(feature_barcodes_ch)

    // add the publish directories to the channel and branch based on existing rad files
    feature_ch = feature_channel
      .map{
        def meta = it.clone();
        meta.feature_rad_publish_dir = "${params.checkpoints_dir}/rad/${meta.library_id}";
        meta.feature_rad_dir = "${meta.feature_rad_publish_dir}/${meta.run_id}-features";
        meta.barcode_file = "${params.barcode_dir}/${cell_barcodes[meta.technology]}";
        meta // return modified meta object
      }
      // branch based on whether mapping should be run (make_rad) or skipped (has_rad)
      // if neither fastq or rad dir are present, run goes into missing_inputs branch
      .branch{
        make_rad: (
          // input files exist
          it.files_directory && file(it.files_directory, type: "dir").exists() && (
            // and repeat has been requested
            params.repeat_mapping
            // or the feature rad file directory does not exist
            || !file(it.feature_rad_dir).exists()
          )
        )
        has_rad: file(it.feature_rad_dir).exists()
        missing_inputs: true
      }

    // send run ids in feature_ch.missing_inputs to log
    feature_ch.missing_inputs
      .subscribe{
        log.error("The expected feature input fastq or rad files for ${it.run_id} are missing.")
      }

    // pull out files that need to be repeated
    feature_reads_ch = feature_ch.make_rad
      // create tuple of [metadata, [Read1 files], [Read2 files]]
      // regex to ensure correct file names if R1 or R2 are in sample identifier
      // We start by including the feature_barcode file so we can combine with the indices, but that will be removed
      .map{meta -> tuple(
        meta.feature_barcode_file,
        meta,
        files("${meta.files_directory}/*_{R1,R1_*}.fastq.gz", checkIfExists: true)
          .findAll{it.name =~ /_R1(_\d+)?.fastq.gz$/},
        files("${meta.files_directory}/*_{R2,R2_*}.fastq.gz", checkIfExists: true)
          .findAll{it.name =~ /_R2(_\d+)?.fastq.gz$/}
      )}
      .combine(index_feature.out, by: 0) // combine by the feature_barcode_file (reused indices, so combine is needed)
      .map{it.drop(1)} // remove the first element (feature_barcode_file)

    // // if the rad directory has been created and repeat_mapping is set to false
    // create tuple of metadata map (read from output) and rad_directory to be used directly as input to alevin-fry quantification
    feature_rad_ch = feature_ch.has_rad
      .map{meta -> tuple(
        Utils.readMeta(file("${meta.feature_rad_dir}/scpca-meta.json")),
        file(meta.feature_rad_dir, type: 'dir', checkIfExists: true)
      )}

    // run Alevin on feature reads
    alevin_feature(feature_reads_ch)

    // combine output from running alevin step with channel containing libraries that skipped creating RAD file
    all_feature_rad_ch = alevin_feature.out.mix(feature_rad_ch)
      .map{it.toList() + [file(it[0].barcode_file)]}


    // quantify feature reads
    fry_quant_feature(all_feature_rad_ch)

  emit: fry_quant_feature.out
  // a tuple of metadata map and the alevin-fry output directory
}
