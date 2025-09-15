// generate unfiltered and filtered RDS files using scpcaTools

// RNA only libraries
process make_unfiltered_sce {
    container params.SCPCATOOLS_SLIM_CONTAINER
    label 'mem_8'
    tag "${meta.library_id}"
    input:
        tuple val(meta), path(alevin_dir), 
              path(mito_file), path(ref_gtf), 
              path(submitter_cell_types_file), path(openscpca_cell_types_file)
        path sample_metafile
    output:
        tuple val(meta), path(unfiltered_rds)
    script:
        unfiltered_rds = "${meta.library_id}_unfiltered.rds"
        """
        generate_unfiltered_sce_alevin.R \
          --alevin_dir ${alevin_dir} \
          --unfiltered_file ${unfiltered_rds} \
          --technology ${meta.technology} \
          --seq_unit ${meta.seq_unit} \
          --library_id "${meta.library_id}" \
          --sample_id "${meta.sample_id}" \
          --project_id "${meta.project_id}" \
          ${meta.assay_ontology_term_id? "--assay_ontology_term_id ${meta.assay_ontology_term_id}" : ""} \
          ${params.spliced_only ? '--spliced_only' : ''}

        format_unfiltered_sce.R \
          --sce_file ${unfiltered_rds} \
          --mito_file ${mito_file} \
          --gtf_file ${ref_gtf} \
          --library_id "${meta.library_id}" \
          --sample_id "${meta.sample_id}" \
          --sample_metadata_file ${sample_metafile}

        # Only run scripts if annotations are available:
        if [[ -f "${submitter_cell_types_file}" ]]; then
          add_submitter_annotations.R \
            --sce_file "${unfiltered_rds}" \
            --library_id "${meta.library_id}" \
            --submitter_cell_types_file "${submitter_cell_types_file}"
        fi

        if [[ -n "${openscpca_cell_types_file}" ]]; then
          add_openscpca_annotations.R \
            --sce_file "${unfiltered_rds}" \
            --library_id "${meta.library_id}" \
            --openscpca_cell_types_file "${openscpca_cell_types_file}"
        fi

        """
    stub:
        unfiltered_rds = "${meta.library_id}_unfiltered.rds"
        """
        touch ${unfiltered_rds}
        """
}

// channels with RNA and feature data
process make_unfiltered_sce_with_feature {
    label 'mem_8'
    tag "${rna_meta.library_id}"
    container params.SCPCATOOLS_SLIM_CONTAINER
    input:
        tuple val(feature_meta), path(feature_alevin_dir),
              val(rna_meta), path(alevin_dir),
              path(mito_file), path(ref_gtf), path(submitter_cell_types_file), path(openscpca_cell_types_file)
        path sample_metafile
    output:
        tuple val(meta), path(unfiltered_rds)
    script:
        // add feature metadata as elements of the main meta object
        meta = rna_meta.clone()
        meta['feature_type'] = feature_meta.technology.split('_')[0]
        meta['feature_meta'] = feature_meta

        // If feature_type is "CITEseq", make it "adt"
        if (meta['feature_type'] == "CITEseq") {
          meta['feature_type'] = "adt"
        }

        unfiltered_rds = "${meta.library_id}_unfiltered.rds"
        """
        generate_unfiltered_sce_alevin.R \
          --alevin_dir ${alevin_dir} \
          --feature_dir ${feature_alevin_dir} \
          --feature_name ${meta.feature_type} \
          --unfiltered_file ${unfiltered_rds} \
          --technology ${meta.technology} \
          --seq_unit ${meta.seq_unit} \
          --library_id "${meta.library_id}" \
          --sample_id "${meta.sample_id}" \
          --project_id "${meta.project_id}" \
          ${meta.assay_ontology_term_id? "--assay_ontology_term_id ${meta.assay_ontology_term_id}" : ""} \
          ${params.spliced_only ? '--spliced_only' : ''}

        format_unfiltered_sce.R \
          --sce_file ${unfiltered_rds} \
          --mito_file ${mito_file} \
          --gtf_file ${ref_gtf} \
          --library_id "${meta.library_id}" \
          --sample_id "${meta.sample_id}" \
          --sample_metadata_file ${sample_metafile}

        # Only run script if annotations are available:
        if [[ -n "${submitter_cell_types_file}" ]]; then
          add_submitter_annotations.R \
            --sce_file "${unfiltered_rds}" \
            --library_id "${meta.library_id}" \
            --submitter_cell_types_file "${submitter_cell_types_file}"
        fi

        if [[ -n "${openscpca_cell_types_file}" ]]; then
          add_openscpca_annotations.R \
            --sce_file "${unfiltered_rds}" \
            --library_id "${meta.library_id}" \
            --openscpca_cell_types_file "${openscpca_cell_types_file}"
        fi
        """
    stub:
        meta = rna_meta.clone()
        meta['feature_type'] = feature_meta.technology.split('_')[0]
        meta['feature_meta'] = feature_meta

        unfiltered_rds = "${meta.library_id}_unfiltered.rds"
        """
        touch "${meta.library_id}_unfiltered.rds"
        """
}

process make_unfiltered_sce_cellranger {
    container params.SCPCATOOLS_SLIM_CONTAINER
    label 'mem_8'
    tag "${meta.library_id}"
    input:
        tuple val(meta), path(cellranger_dir), 
              path(versions_file), path(metrics_file), 
              path(ref_gtf), path(submitter_cell_types_file), path(openscpca_cell_types_file)
        path sample_metafile
    output:
        tuple val(meta), path(unfiltered_rds)
    script:
        unfiltered_rds = "${meta.library_id}_unfiltered.rds"
        """
        generate_unfiltered_sce_cellranger.R \
          --cellranger_dir ${cellranger_dir} \
          --unfiltered_file ${unfiltered_rds} \
          --versions_file ${versions_file} \
          --metrics_file ${metrics_file} \
          --reference_index ${meta.cellranger_index} \
          --reference_probeset ${meta.flex_probeset} \
          --technology ${meta.technology} \
          --seq_unit ${meta.seq_unit} \
          --library_id "${meta.library_id}" \
          --sample_id "${meta.sample_id}" \
          --project_id "${meta.project_id}" \
          ${meta.assay_ontology_term_id? "--assay_ontology_term_id ${meta.assay_ontology_term_id}" : ""}

        format_unfiltered_sce.R \
          --sce_file ${unfiltered_rds} \
          --gtf_file ${ref_gtf} \
          --library_id "${meta.library_id}" \
          --sample_id "${meta.sample_id}" \
          --sample_metadata_file ${sample_metafile}

        # Only run script if annotations are available:
        if [[ -n "${submitter_cell_types_file}" ]]; then
          add_submitter_annotations.R \
            --sce_file "${unfiltered_rds}" \
            --library_id "${meta.library_id}" \
            --submitter_cell_types_file "${submitter_cell_types_file}"
        fi

        if [[ -n "${openscpca_cell_types_file}" ]]; then
          add_openscpca_annotations.R \
            --sce_file "${unfiltered_rds}" \
            --library_id "${meta.library_id}" \
            --openscpca_cell_types_file "${openscpca_cell_types_file}"
        fi

        """
    stub:
        unfiltered_rds = "${meta.library_id}_unfiltered.rds"
        """
        touch ${unfiltered_rds}
        """
}

process filter_sce {
  container params.SCPCATOOLS_SLIM_CONTAINER
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
      --filtered_file "filtered.rds" \
      ${adt_present ? "--adt_name ${meta.feature_type}":""} \
      ${adt_present ? "--adt_barcode_file ${feature_barcode_file}":""} \
      --prob_compromised_cutoff ${params.prob_compromised_cutoff} \
      ${params.seed ? "--random_seed ${params.seed}" : ""} \
      --no_sce_compression

    detect_doublets.R \
      --input_sce_file "filtered.rds" \
      --output_sce_file ${filtered_rds} \
      ${params.seed ? "--random_seed ${params.seed}" : ""} \
      --threads ${task.cpus}
    """
  stub:
    filtered_rds = "${meta.library_id}_filtered.rds"
    """
    touch ${filtered_rds}
    """
}

process genetic_demux_sce {
  container params.SCPCATOOLS_SLIM_CONTAINER
  label 'mem_8'
  tag "${meta.library_id}"
  input:
    tuple val(demux_meta), path(vireo_dir),
          val(meta), path(unfiltered_rds), path(filtered_rds)
  output:
    tuple val(meta), path(unfiltered_rds), path(demux_rds)
  script:
    // demultiplex results will be added to the SCE object colData
    demux_rds = "${filtered_rds.baseName}_genetic-demux.rds"
    """
    add_demux_sce.R \
      --sce_file ${filtered_rds} \
      --output_sce_file ${demux_rds} \
      --library_id ${meta.library_id} \
      --vireo_dir ${vireo_dir}
    """
  stub:
    demux_rds = "${filtered_rds.baseName}_genetic-demux.rds"
    """
    touch ${demux_rds}
    """
}

process cellhash_demux_sce {
  container params.SCPCATOOLS_SEURAT_CONTAINER
  label 'mem_8'
  tag "${meta.library_id}"
  input:
    tuple val(meta), path(unfiltered_rds), path(filtered_rds)
    path cellhash_pool_file
  output:
    tuple val(meta), path(unfiltered_rds), path(demux_rds)
  script:
    // demultiplex results will be added to the SCE object colData
    demux_rds = "${filtered_rds.baseName}_cellhash-demux.rds"
    """
    add_demux_sce.R \
      --sce_file ${filtered_rds}  \
      --output_sce_file ${demux_rds} \
      --library_id ${meta.library_id} \
      --cellhash_pool_file ${cellhash_pool_file} \
      --hash_demux \
      --seurat_demux
    """
  stub:
    demux_rds = "${filtered_rds.baseName}_cellhash-demux.rds"
    """
    touch ${demux_rds}
    """
}

process post_process_sce {
  container params.SCPCATOOLS_SLIM_CONTAINER
  label 'mem_8'
  tag "${meta.library_id}"
  input:
    tuple val(meta), path(unfiltered_rds), path(filtered_rds)
  output:
    tuple val(meta), path(unfiltered_rds), path(filter_labeled_rds), path(processed_rds)
  script:
    filter_labeled_rds = "${meta.library_id}_filtered_labeled.rds"
    processed_rds = "${meta.library_id}_processed.rds"
    """
    post_process_sce.R \
      --filtered_sce_file ${filtered_rds} \
      --out_filtered_sce_file ${filter_labeled_rds} \
      --out_processed_sce_file ${processed_rds} \
      --gene_cutoff ${params.gene_cutoff} \
      --n_hvg ${params.num_hvg} \
      --n_pcs ${params.num_pcs} \
      ${params.seed ? "--random_seed ${params.seed}" : ""}
    """
  stub:
    filter_labeled_rds = "${meta.library_id}_filtered_labeled.rds"
    processed_rds = "${meta.library_id}_processed.rds"
    """
    touch ${filter_labeled_rds}
    touch ${processed_rds}
    """
}




workflow generate_sce {
  // generate rds files for RNA-only samples
  take:
    quant_channel
    sample_metafile
  main:
    def empty_file = "${projectDir}/assets/NO_FILE"
    
    sce_ch = quant_channel
      .map{it.toList() + [file(it[0].mito_file, checkIfExists: true),
                          file(it[0].ref_gtf, checkIfExists: true),
                          // either submitter/openscpca cell type files, or empty array if not available
                          it[0].submitter_cell_types_file ? file(it[0].submitter_cell_types_file, checkIfExists: true) : [],
                          it[0].openscpca_cell_types_file ? file(it[0].openscpca_cell_types_file, checkIfExists: true) : []
                         ]}

    make_unfiltered_sce(sce_ch, sample_metafile)

    // provide empty feature barcode file, since no features here
    unfiltered_sce_ch = make_unfiltered_sce.out
      .map{it.toList() + [file(empty_file, checkIfExists: true)]}

    filter_sce(unfiltered_sce_ch)

  emit: filter_sce.out
  // a tuple of meta and the filtered and unfiltered rds files
}

workflow generate_sce_with_feature {
  // generate rds files for feature + quant samples
  // input is a channel with feature_meta, feature_quantdir, rna_meta, rna_quantdir
  take:
    feature_quant_channel
    sample_metafile
  main:
    def empty_file = "${projectDir}/assets/NO_FILE"

    feature_sce_ch = feature_quant_channel
      // RNA meta is in the third slot here
      .map{it.toList() + [file(it[2].mito_file, checkIfExists: true),
                          file(it[2].ref_gtf, checkIfExists: true),
                          // either submitter/openscpca cell type files, or empty array if not available
                          it[0].submitter_cell_types_file ? file(it[0].submitter_cell_types_file, checkIfExists: true) : [],
                          it[0].openscpca_cell_types_file ? file(it[0].openscpca_cell_types_file, checkIfExists: true) : []
                         ]}

    make_unfiltered_sce_with_feature(feature_sce_ch, sample_metafile)

    // append the feature barcode file
    unfiltered_feature_sce_ch = make_unfiltered_sce_with_feature.out
      .map{it.toList() + [file(it[0]["feature_meta"].feature_barcode_file ?: empty_file, checkIfExists: true)]}

    filter_sce(unfiltered_feature_sce_ch)

  emit: filter_sce.out
  // a tuple of meta and the filtered and unfiltered rds files
}

workflow generate_sce_cellranger {
  // generate rds files for RNA-only samples from cellranger multi
  take:
    quant_channel
    sample_metafile
  main:
    def empty_file = "${projectDir}/assets/NO_FILE"

    sce_ch = quant_channel
      .map{it.toList() + [file(it[0].ref_gtf, checkIfExists: true),
                          // either submitter/openscpca cell type files, or empty array if not available
                          it[0].submitter_cell_types_file ? file(it[0].submitter_cell_types_file, checkIfExists: true) : [],
                          it[0].openscpca_cell_types_file ? file(it[0].openscpca_cell_types_file, checkIfExists: true) : []
                         ]}

    make_unfiltered_sce_cellranger(sce_ch, sample_metafile)

    // provide empty feature barcode file, since no features here
    unfiltered_sce_ch = make_unfiltered_sce_cellranger.out
      .map{it.toList() + [file(empty_file, checkIfExists: true)]}

    filter_sce(unfiltered_sce_ch)

  emit: filter_sce.out
  // a tuple of meta and the filtered and unfiltered rds files
}
