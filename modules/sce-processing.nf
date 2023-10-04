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
          --seq_unit ${meta.seq_unit} \
          --library_id "${meta.library_id}" \
          --sample_id "${meta.sample_id}" \
          --project_id "${meta.project_id}" \
          --sample_metadata_file ${sample_metafile} \
          ${meta.assay_ontology_term_id? "--assay_ontology_term_id ${meta.assay_ontology_term_id}" : ""} \
          ${params.spliced_only ? '--spliced_only' : ''}


        # Only run script if annotations are available:
        if [ "${submitter_cell_types_file.name}" != "NO_FILE" ]; then
          add_submitter_annotations.R \
            --sce_file "${unfiltered_rds}" \
            --library_id "${meta.library_id}" \
            --submitter_cell_types_file "${submitter_cell_types_file}"
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
              val (rna_meta), path(alevin_dir),
              path(mito_file), path(ref_gtf), path(submitter_cell_types_file)
        path sample_metafile
    output:
        tuple val(meta), path(unfiltered_rds)
    script:
        // add feature metadata as elements of the main meta object
        meta = rna_meta
        meta['feature_type'] = feature_meta.technology.split('_')[0]
        meta['feature_technology'] = feature_meta.technology
        meta['feature_run_id'] = feature_meta.run_id
        meta['feature_library_id'] = feature_meta.library_id
        meta['feature_sample_id'] = feature_meta.sample_id
        meta['feature_project_id'] = feature_meta.project_id
        meta['feature_barcode_file'] = feature_meta.feature_barcode_file
        meta['feature_barcode_geom'] = feature_meta.feature_barcode_geom
        meta['feature_files_directory'] = feature_meta.files_directory

        // If feature_type is "CITEseq", make it "adt"
        if (meta['feature_type'] == "CITEseq") {
          meta['feature_type'] = "adt"
        }

        unfiltered_rds = "${meta.library_id}_unfiltered.rds"
        """
        generate_unfiltered_sce.R \
          --alevin_dir ${alevin_dir} \
          --feature_dir ${feature_alevin_dir} \
          --feature_name ${meta.feature_type} \
          --unfiltered_file ${unfiltered_rds} \
          --mito_file ${mito_file} \
          --gtf_file ${ref_gtf} \
          --technology ${meta.technology} \
          --seq_unit ${meta.seq_unit} \
          --library_id "${meta.library_id}" \
          --sample_id "${meta.sample_id}" \
          --project_id "${meta.project_id}" \
          --sample_metadata_file ${sample_metafile} \
          ${meta.assay_ontology_term_id? "--assay_ontology_term_id ${meta.assay_ontology_term_id}" : ""} \
          ${params.spliced_only ? '--spliced_only' : ''}

        # Only run script if annotations are available:
        if [ ${submitter_cell_types_file.name} != "NO_FILE" ]; then
          add_submitter_annotations.R \
            --sce_file "${unfiltered_rds}" \
            --library_id "${meta.library_id}" \
            --submitter_cell_types_file "${submitter_cell_types_file}"
        fi

        """
    stub:
        meta = rna_meta
        meta['feature_type'] = feature_meta.technology.split('_')[0]
        meta['feature_technology'] = feature_meta.technology
        meta['feature_run_id'] = feature_meta.run_id
        meta['feature_library_id'] = feature_meta.library_id
        meta['feature_sample_id'] = feature_meta.sample_id
        meta['feature_project_id'] = feature_meta.project_id
        meta['feature_barcode_file'] = feature_meta.feature_barcode_file ?: ""
        meta['feature_barcode_geom'] = feature_meta.feature_barcode_geom ?: ""
        meta['feature_files_directory'] = feature_meta.files_directory

        unfiltered_rds = "${meta.library_id}_unfiltered.rds"
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

        // Checks for whether we have ADT data:
        // - feature_type should be adt
        // - barcode file should _not_ be the empty file NO_FILE
        adt_present = meta.feature_type == 'adt' &
          feature_barcode_file.name != "NO_FILE"

        """
        filter_sce.R \
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
    tuple val(meta), path(unfiltered_rds), path(filtered_demux_rds)
  script:
    // output will be same as input, with replacement of the filtered_rds file
    // demultiplex results will be added to the SCE object colData
    filtered_demux_rds = "${meta.library_id}_demux_filtered.rds"
    """
    add_demux_sce.R \
      --sce_file ${filtered_rds} \
      --output_sce_file ${filtered_demux_rds} \
      --library_id ${meta.library_id} \
      --vireo_dir ${vireo_dir}
    """
  stub:
    filtered_demux_rds = "${meta.library_id}_demux_filtered.rds"
    """
    touch ${filtered_demux_rds}
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
    tuple val(meta), path(unfiltered_rds), path(demux_rds)
  script:
    // output will be same as input, with replacement of the filtered_rds file
    // demultiplex results will be added to the SCE object colData
    demux_rds = "${meta.library_id}_demux_filtered.rds"
    """
    add_demux_sce.R \
      --sce_file "${filtered_rds}" \
      --output_sce_file "${demux_rds}" \
      --library_id ${meta.library_id} \
      --cellhash_pool_file "${cellhash_pool_file}" \
      --hash_demux \
      --seurat_demux
    """
  stub:
    demux_rds = "${meta.library_id}_demux_filtered.rds"
    """
    touch ${demux_rds}
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


// used when a given file is not defined in the below workflows
empty_file = "${projectDir}/assets/NO_FILE"

workflow generate_sce {
  // generate rds files for RNA-only samples
  take:
    quant_channel
    sample_metafile
  main:

    sce_ch = quant_channel
      .map{it.toList() + [file(it[0].mito_file),
                          file(it[0].ref_gtf),
                          // either submitter cell type files, or empty file if not available
                          file(it[0].submitter_cell_types_file ?: empty_file)
                         ]}

    make_unfiltered_sce(sce_ch, sample_metafile)

    // provide empty feature barcode file, since no features here
    unfiltered_sce_ch = make_unfiltered_sce.out
      .map{it.toList() + [file(empty_file)]}

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

    feature_sce_ch = feature_quant_channel
      // RNA meta is in the third slot here
      .map{it.toList() + [file(it[2].mito_file),
                          file(it[2].ref_gtf),
                          // either submitter cell type files, or empty file if not available
                          file(it[2].submitter_cell_types_file ?: empty_file)
                         ]}

    make_merged_unfiltered_sce(feature_sce_ch, sample_metafile)

    // append the feature barcode file
    unfiltered_merged_sce_ch = make_merged_unfiltered_sce.out
      .map{it.toList() + [file(it[0]["feature_barcode_file"] ?: empty_file)]}

    filter_sce(unfiltered_merged_sce_ch)

  emit: filter_sce.out
  // a tuple of meta and the filtered and unfiltered rds files
}
