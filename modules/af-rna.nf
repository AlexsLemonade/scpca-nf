

// generates RAD file using alevin
process alevin_rad{
  container params.SALMON_CONTAINER
  label 'cpus_12'
  label 'mem_24'
  label 'disk_big'
  tag "${meta.run_id}-rna"
  publishDir "${meta.rad_publish_dir}"
  input:
    tuple val(meta),
          path(read1), path(read2)
    path index
  output:
    tuple val(meta), path(rad_dir)
  script:
    // define the local location for the rad output
    rad_dir = file(meta.rad_dir).name
    // choose flag by technology
    tech_flag = ['10Xv2': '--chromium',
                 '10Xv2_5prime': '--chromium',
                 '10Xv3': '--chromiumV3',
                 '10Xv3.1': '--chromiumV3']
    // get meta to write as file
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
}

// quantify rna from RAD input
process fry_quant_rna{
  container params.ALEVINFRY_CONTAINER
  label 'cpus_8'
  label 'mem_8'
  tag "${meta.run_id}-rna"
  publishDir "${params.checkpoints_dir}/alevinfry/${meta.library_id}", enabled: params.publish_fry_outs

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
    // add rad publish directory, rad directory, and barcode file to meta
    rna_channel = rna_channel
      .map{it.rad_publish_dir = "${params.checkpoints_dir}/rad/${it.library_id}";
           it.rad_dir = "${it.rad_publish_dir}/${it.run_id}-rna";
           it.barcode_file = "${params.barcode_dir}/${params.cell_barcodes[it.technology]}";
           it}
       // split based in whether repeat_mapping is false and a previous dir exists
      .branch{
          has_rad: !params.repeat_mapping && file(it.rad_dir).exists()
          make_rad: true
       }

    // If We need to create rad files, create a new channel with tuple of (metadata map, [Read1 files], [Read2 files])
    rna_reads_ch = rna_channel.make_rad
      .map{meta -> tuple(meta,
                         file("${meta.files_directory}/*_R1_*.fastq.gz"),
                         file("${meta.files_directory}/*_R2_*.fastq.gz")
                        )}

    // if the rad directory has been created and repeat_mapping is set to false
    // create tuple of (metdata map, rad_directory) to be used directly as input to alevin-fry quantification
    rna_rad_ch = rna_channel.has_rad
      .map{meta -> tuple(Utils.readMeta(file("${meta.rad_dir}/scpca-meta.json")),
                         file(meta.rad_dir)
                         )}

    // run Alevin for mapping on libraries that don't have RAD directory already created
    alevin_rad(rna_reads_ch, params.splici_index)

    // combine ouput from running alevin step with channel containing libraries that skipped creating a RAD file
    all_rad_ch = alevin_rad.out.mix(rna_rad_ch)
      .map{it.toList() << file(it[0].barcode_file)} // add barcode file to channel to use in fry_quant_rna process

    // quantify with alevin-fry
    fry_quant_rna(all_rad_ch, params.t2g_3col_path)

  emit: fry_quant_rna.out
  // a tuple of meta and the alevin-fry output directory
}
