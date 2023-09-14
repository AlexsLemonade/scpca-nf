// generate unfiltered and filtered RDS files using scpcaTools

// RNA only libraries
process make_unfiltered_sce{
    container params.SCPCATOOLS_CONTAINER
    label 'mem_8'
    tag "${meta.library_id}"
    input:
        tuple val(meta), path(alevin_dir), path(mito_file), path(ref_gtf), path(submitter_cell_types_file)
        path sample_metafile
    output:
        tuple val(meta), path(unfiltered_rds)
    script:
        unfiltered_rds = "${meta.library_id}_unfiltered.rds"
      
        """
        generate_unfiltered_sce.R \
          --alevin_dir ${alevin_dir} \
          --unfiltered_file ${unfiltered_rds} \
          --mito_file ${mito_file} \
          --gtf_file ${ref_gtf} \
          --technology ${meta.technology} \
          --library_id "${meta.library_id}" \
          --sample_id "${meta.sample_id}" \
          --sample_metadata_file ${sample_metafile} \
          ${params.spliced_only ? '--spliced_only' : ''}
          
          
        # Only run script if annotations are available:
        if [ ${submitter_cell_types_file.name} != "NO_FILE.txt" ]; then
          add_submitter_annotations.R \
            --unfiltered_file ${unfiltered_rds} \
            --library_id "${meta.library_id}" \
            --submitter_cell_types_file ${submitter_cell_types_file}
        fi

        """
    stub:
        unfiltered_rds = "${meta.library_id}_unfiltered.rds"
        """
        touch ${unfiltered_rds}
        """
}

// channels with RNA and feature data
process make_merged_unfiltered_sce{
    label 'mem_8'
    tag "${meta.library_id}"
    container params.SCPCATOOLS_CONTAINER
    input:
        tuple val(feature_meta), path(feature_alevin_dir),
              val (meta), path(alevin_dir),
              path(mito_file), path(ref_gtf), path(submitter_cell_types_file)
        path sample_metafile
    output:
        tuple val(meta), path(unfiltered_rds)
    script:
        unfiltered_rds = "${meta.library_id}_unfiltered.rds"
        // add feature metadata as an element of the main meta object
        meta['feature_type'] = feature_meta.technology.split('_')[0]
        meta['feature_meta'] = feature_meta

        // If feature_type is "CITEseq", make it "adt"
        if (meta['feature_type'] == "CITEseq") {
          meta['feature_type'] = "adt"
        }

        """
        generate_unfiltered_sce.R \
          --alevin_dir ${alevin_dir} \
          --feature_dir ${feature_alevin_dir} \
          --feature_name ${meta.feature_type} \
          --unfiltered_file ${unfiltered_rds} \
          --mito_file ${mito_file} \
          --gtf_file ${ref_gtf} \
          --technology ${meta.technology} \
          --library_id "${meta.library_id}" \
          --sample_id "${meta.sample_id}" \
          --sample_metadata_file ${sample_metafile} \
          ${params.spliced_only ? '--spliced_only' : ''}
          
        # Only run script if annotations are available:
        if [ ${submitter_cell_types_file.name} != "NO_FILE.txt" ]; then
          add_submitter_annotations.R \
            --unfiltered_file ${unfiltered_rds} \
            --library_id "${meta.library_id}" \
            --submitter_cell_types_file ${submitter_cell_types_file}
        fi

        """
    stub:
        unfiltered_rds = "${meta.library_id}_unfiltered.rds"
        meta['feature_type'] = feature_meta.technology.split('_')[0]
        meta['feature_meta'] = feature_meta
        """
        touch "${meta.library_id}_unfiltered.rds"
        """
}

process filter_sce{
    container params.SCPCATOOLS_CONTAINER
    label 'mem_8'
    tag "${meta.library_id}"
    input:
        tuple val(meta), path(unfiltered_rds), path(feature_barcode_file)
    output:
        tuple val(meta), path(unfiltered_rds), path(filtered_rds)
    script:
        filtered_rds = "${meta.library_id}_filtered.rds"

        // Three checks for whether we have ADT data:
        // - feature_type should be adt
        // - barcode file should exist
        // - barcode file should _not_ be the empty file NO_FILE.txt
        adt_present = meta.feature_type == 'adt' &
          feature_barcode_file.exists() &
          feature_barcode_file.name != "NO_FILE.txt"

        """
        filter_sce_rds.R \
          --unfiltered_file ${unfiltered_rds} \
          --filtered_file ${filtered_rds} \
          ${adt_present ? "--adt_name ${meta.feature_type}":""} \
          ${adt_present ? "--adt_barcode_file ${feature_barcode_file}":""} \
          --prob_compromised_cutoff ${params.prob_compromised_cutoff} \
          ${params.seed ? "--random_seed ${params.seed}" : ""}
        """
    stub:
        filtered_rds = "${meta.library_id}_filtered.rds"
        """
        touch ${filtered_rds}
        """
}

process genetic_demux_sce{
  container params.SCPCATOOLS_CONTAINER
  label 'mem_8'
  tag "${meta.library_id}"
  input:
    tuple val(demux_meta), path(vireo_dir),
          val(meta), path(unfiltered_rds), path(filtered_rds)
  output:
    tuple val(meta), path(unfiltered_rds), path(filtered_rds)
  script:
    // output will be same as input, with replacement of the filtered_rds file
    // demultiplex results will be added to the SCE object colData
    """
    mv ${filtered_rds} filtered_nodemux.rds
    add_demux_sce.R \
      --sce_file filtered_nodemux.rds \
      --output_sce_file ${filtered_rds} \
      --library_id ${meta.library_id} \
      --vireo_dir ${vireo_dir}
    """
  stub:
    """
    touch ${filtered_rds}
    """
}

process cellhash_demux_sce{
  container params.SCPCATOOLS_CONTAINER
  label 'mem_8'
  tag "${meta.library_id}"
  input:
    tuple val(meta), path(unfiltered_rds), path(filtered_rds)
    path cellhash_pool_file
  output:
    tuple val(meta), path(unfiltered_rds), path(filtered_rds)
  script:
    // output will be same as input, with replacement of the filtered_rds file
    // demultiplex results will be added to the SCE object colData
    """
    mv ${filtered_rds} filtered_nodemux.rds
    add_demux_sce.R \
      --sce_file filtered_nodemux.rds \
      --output_sce_file ${filtered_rds} \
      --library_id ${meta.library_id} \
      --cellhash_pool_file ${cellhash_pool_file} \
      --hash_demux \
      --seurat_demux
    """
  stub:
    """
    touch ${filtered_rds}
    """
}

process post_process_sce{
    container params.SCPCATOOLS_CONTAINER
    label 'mem_8'
    tag "${meta.library_id}"
    input:
        tuple val(meta), path(unfiltered_rds), path(filtered_rds)
    output:
        tuple val(meta), path(unfiltered_rds), path(filtered_rds), path(processed_rds)
    script:
        processed_rds = "${meta.library_id}_processed.rds"

        """
        post_process_sce.R \
          --filtered_sce_file ${filtered_rds} \
          --output_sce_file ${processed_rds} \
          --gene_cutoff ${params.gene_cutoff} \
          --n_hvg ${params.num_hvg} \
          --n_pcs ${params.num_pcs} \
          ${params.seed ? "--random_seed ${params.seed}" : ""}
        """
    stub:
        processed_rds = "${meta.library_id}_processed.rds"
        """
        touch ${processed_rds}
        """
}

workflow generate_sce {
  // generate rds files for RNA-only samples
  take:
    quant_channel
    sample_metafile
  main:
    
    empty_file = file("${projectDir}/assets/NO_FILE.txt")
 
    sce_ch = quant_channel
      .map{it.toList() + [file(it[0].mito_file), 
                         file(it[0].ref_gtf), 
                         // either submitter cell type files, or empty file if not available
                         file(it[0].submitter_cell_types_file).exists() ? file(it[0].submitter_cell_types_file) : empty_file
                        ]}
                 
    make_unfiltered_sce(sce_ch, sample_metafile)
 
    // provide empty feature barcode file, since no features here
    unfiltered_sce_ch = make_unfiltered_sce.out
      .map{it.toList() + [empty_file]}

    filter_sce(unfiltered_sce_ch)

  emit: filter_sce.out
  // a tuple of meta and the filtered and unfiltered rds files
}

workflow generate_merged_sce {
  // generate rds files for feature + quant samples
  // input is a channel with feature_meta, feature_quantdir, rna_meta, rna_quantdir
  take:
    feature_quant_channel
    sample_metafile
  main:
    empty_file = file("${projectDir}/assets/NO_FILE.txt")
 
    feature_sce_ch = feature_quant_channel 
      // RNA meta is in the third slot here
      .map{it.toList() + [file(it[2].mito_file), 
                         file(it[2].ref_gtf), 
                         // either submitter cell type files, or empty file if not available
                         file(it[2].submitter_cell_types_file).exists() ? file(it[2].submitter_cell_types_file) : empty_file
                        ]}
                 
    make_merged_unfiltered_sce(feature_sce_ch, sample_metafile)

    // append the feature barcode file
    unfiltered_merged_sce_ch = make_merged_unfiltered_sce.out
      .map{it.toList() + [file(it[0]["feature_meta"].feature_barcode_file)]}

    filter_sce(unfiltered_merged_sce_ch)

  emit: filter_sce.out
  // a tuple of meta and the filtered and unfiltered rds files
}
