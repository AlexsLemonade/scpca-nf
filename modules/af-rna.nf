
// generates RAD file using alevin
process alevin_rad{
  container params.SALMON_CONTAINER
  label 'cpus_12'
  label 'mem_24'
  label 'disk_big'
  tag "${meta.run_id}-rna"
  publishDir "${meta.rad_publish_dir}", mode: 'copy'
  input:
    tuple val(meta),
          path(read1), path(read2), path(index)
  output:
    tuple val(meta), path(rad_dir)
  script:
    // define the local location for the rad output
    rad_dir = file(meta.rad_dir).name
    // choose flag by technology
    tech_flag = [
      '10Xv2': '--chromium',
      '10Xv2_5prime': '--chromium',
      '10Xv3': '--chromiumV3',
      '10Xv3.1': '--chromiumV3'
    ]
    // get meta to write as file
    meta += Utils.getVersions(workflow, nextflow)
    meta_json = Utils.makeJson(meta)
    // run alevin like normal with the --rad flag
    // creates output directory with RAD file needed for alevin-fry
    """
    mkdir -p ${rad_dir}
    salmon alevin \
      -l ISR \
      ${tech_flag[meta.technology]} \
      -1 ${read1} \
      -2 ${read2} \
      -i ${index} \
      -o ${rad_dir} \
      -p ${task.cpus} \
      --dumpFeatures \
      --rad

    echo '${meta_json}' > ${rad_dir}/scpca-meta.json
    """
  stub:
    rad_dir = file(meta.rad_dir).name
    meta += Utils.getVersions(workflow, nextflow)
    meta_json = Utils.makeJson(meta)
    """
    mkdir -p ${rad_dir}
    echo '${meta_json}' > ${rad_dir}/scpca-meta.json
    """
}

// quantify rna from RAD input
process fry_quant_rna{
  container params.ALEVINFRY_CONTAINER
  label 'cpus_8'
  label 'mem_8'
  tag "${meta.run_id}-rna"
  publishDir "${params.checkpoints_dir}/alevinfry/${meta.library_id}", mode: 'copy', enabled: params.publish_fry_outs

  input:
    tuple val(meta), path(rad_dir), path(barcode_file), path(t2g_3col)
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
      --expected-ori ${meta.technology == '10Xv2_5prime' ? 'rc' : 'fw'} \
      -o ${quant_dir} \
      --unfiltered-pl ${barcode_file}

    alevin-fry collate \
      --input-dir ${quant_dir} \
      --rad-dir ${rad_dir} \
      -t ${task.cpus}

    alevin-fry quant \
      --input-dir ${quant_dir} \
      --tg-map ${t2g_3col} \
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


workflow map_quant_rna {
  take: rna_channel
  // a channel with a map of metadata for each rna library to process
  main:
    // add rad publish directory, rad directory, and barcode file to meta
    rna_channel = rna_channel
      .map{
        def meta = it.clone();
        meta.rad_publish_dir = "${params.checkpoints_dir}/rad/${meta.library_id}";
        meta.rad_dir = "${meta.rad_publish_dir}/${meta.run_id}-rna";
        meta.barcode_file = "${params.barcode_dir}/${params.cell_barcodes[meta.technology]}";
        meta // return modified meta object
      }
       // branch based on whether mapping should be run (make_rad) or skipped (has_rad)
      .branch{
        make_rad: (
          // repeat has been requested
          params.repeat_mapping
          // the rad directory does not exist
          || !file(it.rad_dir).exists()
          // the assembly has changed; if rad_dir doesn't exist, this line won't get hit
          || Utils.getMetaVal(file("${it.rad_dir}/scpca-meta.json"), "ref_assembly") != "${it.ref_assembly}"
        )
        has_rad: true
      }

    // If we need to create rad files, create a new channel with tuple of (metadata map, [Read1 files], [Read2 files])
    rna_reads_ch = rna_channel.make_rad
      .map{meta -> tuple(
        meta,
        // fail if the fastq files do not exist
        file("${meta.files_directory}/*_{R1,R1_*}.fastq.gz", checkIfExists: true),
        file("${meta.files_directory}/*_{R2,R2_*}.fastq.gz", checkIfExists: true),
        file(meta.salmon_splici_index, type: 'dir')
      )}

    // if the rad directory has been created and repeat_mapping is set to false
    // create tuple of metdata map (read from output) and rad_directory to be used directly as input to alevin-fry quantification
    rna_rad_ch = rna_channel.has_rad
      .map{meta -> tuple(
        Utils.readMeta(file("${meta.rad_dir}/scpca-meta.json", checkIfExists: true)),
        file(meta.rad_dir, type: 'dir', checkIfExists: true) // fail if no rad directory
      )}

    // run Alevin for mapping on libraries that don't have RAD directory already created
    alevin_rad(rna_reads_ch)

    // combine output from running alevin step with channel containing libraries that skipped creating a RAD file
    all_rad_ch = alevin_rad.out.mix(rna_rad_ch)
      // add barcode and t2g files to channel to use in fry_quant_rna process
      .map{it.toList() + [file(it[0].barcode_file), file(it[0].t2g_3col_path)]}

    // quantify with alevin-fry
    fry_quant_rna(all_rad_ch)

  emit: fry_quant_rna.out
  // a tuple of meta and the alevin-fry output directory
}
