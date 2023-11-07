
process classify_singler {
  container params.SCPCATOOLS_CONTAINER
  publishDir (
    path: "${meta.celltype_publish_dir}",
    mode: 'copy',
    pattern: "${singler_dir}"
  )
  label 'mem_8'
  label 'cpus_4'
  tag "${meta.library_id}"
  input:
    tuple val(meta), path(processed_rds), path(singler_model_file)
  output:
    tuple val(meta.library_id), path(singler_dir)
  script:
    singler_dir = file(meta.singler_dir).name
    """
    # create output directory
    mkdir "${singler_dir}"

    classify_SingleR.R \
      --sce_file "${processed_rds}" \
      --singler_model_file "${singler_model_file}" \
      --output_singler_results_file "${singler_dir}/singler_results.rds" \
      --seed ${params.seed} \
      --threads ${task.cpus}

    # write out meta file
    echo '${Utils.makeJson(meta)}' > "${singler_dir}/scpca-meta.json"
    """
  stub:
    singler_dir = file(meta.singler_dir).name
    """
    mkdir "${singler_dir}"
    echo '${Utils.makeJson(meta)}' > "${singler_dir}/scpca-meta.json"
    """
}


process classify_cellassign {
  container params.SCPCATOOLS_CONTAINER
    publishDir (
      path: "${meta.celltype_publish_dir}",
      mode: 'copy',
      pattern: "${cellassign_dir}"
    )
  label 'mem_32'
  label 'cpus_12'
  tag "${meta.library_id}"
  input:
    tuple val(meta), path(processed_rds), path(cellassign_reference_file)
  output:
    tuple val(meta.library_id), path(cellassign_dir)
  script:
    cellassign_dir = file(meta.cellassign_dir).name

    """
    # create output directory
    mkdir "${cellassign_dir}"

    # Convert SCE to AnnData
    sce_to_anndata.R \
        --input_sce_file "${processed_rds}" \
        --output_rna_h5 "processed.hdf5"

    # Run CellAssign
    predict_cellassign.py \
      --input_hdf5_file "processed.hdf5" \
      --output_predictions "${cellassign_dir}/cellassign_predictions.tsv" \
      --reference "${cellassign_reference_file}" \
      --seed ${params.seed} \
      --threads ${task.cpus}

    # write out meta file
    echo '${Utils.makeJson(meta)}' > "${cellassign_dir}/scpca-meta.json"
    """
  stub:
    cellassign_dir = file(meta.cellassign_dir).name
    """
    mkdir "${cellassign_dir}"
    echo '${Utils.makeJson(meta)}' > "${cellassign_dir}/scpca-meta.json"
    """
}

process add_celltypes_to_sce {
  container params.SCPCATOOLS_CONTAINER
  label 'mem_4'
  label 'cpus_2'
  tag "${meta.library_id}"
  input:
    tuple val(meta), path(processed_rds), path(singler_dir), path(cellassign_dir)
  output:
    tuple val(meta), path(annotated_rds)
  script:
    annotated_rds = "${meta.library_id}_processed_annotated.rds"
    singler_present = "${singler_dir.name}" != "NO_FILE"
    singler_results = "${singler_dir}/singler_results.rds"
    cellassign_present = "${cellassign_dir.name}" != "NO_FILE"
    cellassign_predictions = "${cellassign_dir}/cellassign_predictions.tsv"
    cellassign_ref_name = file("${meta.cellassign_reference_file}").baseName
    """
    add_celltypes_to_sce.R \
      --input_sce_file ${processed_rds} \
      --output_sce_file ${annotated_rds} \
      ${singler_present ? "--singler_results  ${singler_results}" : ''} \
      ${cellassign_present ? "--cellassign_predictions  ${cellassign_predictions}" : ''} \
      ${cellassign_present ? "--cellassign_ref_name ${cellassign_ref_name}" : ''}
    """
  stub:
    annotated_rds = "${meta.library_id}_processed_annotated.rds"
    """
    touch ${annotated_rds}
    """
}

empty_file = "${projectDir}/assets/NO_FILE"

workflow annotate_celltypes {
  take: sce_files_channel // channel of meta, unfiltered_sce, filtered_sce, processed_sce
  main:
    // get just the meta and processed sce
    processed_sce_channel = sce_files_channel.map{[it[0], it[3]]}

    // channel with celltype model and project ids
    celltype_ch = Channel.fromPath(params.celltype_project_metafile)
      .splitCsv(header: true, sep: '\t')
      .map{[
        // project id
        it.scpca_project_id,
        // singler model file
        Utils.parseNA(it.singler_ref_file) ? "${params.singler_models_dir}/${it.singler_ref_file}" : '',
        // cellassign reference file
        Utils.parseNA(it.cellassign_ref_file) ? "${params.cellassign_ref_dir}/${it.cellassign_ref_file}" : ''
      ]}

    // create input for typing: [augmented meta, processed_sce]
    celltype_input_ch = processed_sce_channel
      .map{[it[0]["project_id"]] + it}
      .combine(celltype_ch, by: 0)
      // current contents: [project_id, meta, processed_sce, singler_model_file, cellassign_reference_file]
      // add values to meta for later use
      .map{ project_id, meta_in, processed_sce, singler_model_file, cellassign_reference_file ->
        def meta = meta_in.clone(); // local copy for safe modification
        meta.celltype_publish_dir = "${params.checkpoints_dir}/celltype/${meta.library_id}";
        meta.singler_dir = "${meta.celltype_publish_dir}/${meta.library_id}_singler";
        meta.cellassign_dir = "${meta.celltype_publish_dir}/${meta.library_id}_cellassign";
        meta.singler_model_file = singler_model_file;
        meta.cellassign_reference_file = cellassign_reference_file;
        meta.singler_results_file = "${meta.singler_dir}/singler_results.rds";
        meta.cellassign_predictions_file = "${meta.cellassign_dir}/cellassign_predictions.tsv"
        // return simplified input:
        [meta, processed_sce]
      }


    // creates [meta, processed sce, singler model file]
    singler_input_ch = celltype_input_ch
      // add in singler model or empty file
      .map{it.toList() + [file(it[0].singler_model_file ?: empty_file)]}
      // skip if no singleR model file or if singleR results are already present
      .branch{
        skip_singler: !params.repeat_celltyping && file(it[0].singler_results_file).exists()
        missing_ref: it[2].name == "NO_FILE"
        do_singler: true
      }


    // perform singleR celltyping and export results
    classify_singler(singler_input_ch.do_singler)

    // singleR output channel: [library_id, singler_results]
    singler_output_ch = singler_input_ch.skip_singler
      // provide existing singler results dir for those we skipped
      .map{[it[0]["library_id"], file(it[0].singler_dir, type: 'dir')]}
      // add empty file for missing ref samples
      .mix(singler_input_ch.missing_ref.map{[it[0]["library_id"], file(empty_file)]} )
      // add in channel outputs
      .mix(classify_singler.out)

    // create cellassign input channel: [meta, processed sce, cellassign reference file]
      cellassign_input_ch = celltype_input_ch
      // add in cellassign reference
      .map{it.toList() + [file(it[0].cellassign_reference_file ?: empty_file)]}
      // skip if no cellassign reference file or reference name is not defined
      .branch{
        skip_cellassign: !params.repeat_celltyping && file(it[0].cellassign_predictions_file).exists()
        missing_ref: it[2].name == "NO_FILE"
        do_cellassign: true
      }


    // perform CellAssign celltyping and export results
    classify_cellassign(cellassign_input_ch.do_cellassign)

    // cellassign output channel: [library_id, cellassign_dir]
    cellassign_output_ch = cellassign_input_ch.skip_cellassign
      // provide existing cellassign predictions dir for those we skipped
      .map{[it[0]["library_id"], file(it[0].cellassign_dir, type: 'dir')]}
      // add empty file for missing ref samples
      .mix(cellassign_input_ch.missing_ref.map{[it[0]["library_id"], file(empty_file)]} )
      // add in channel outputs
      .mix(classify_cellassign.out)

    // prepare input for process to add celltypes to the processed SCE
    // result is [meta, processed rds, singler dir, cellassign dir]
    assignment_input_ch = celltype_input_ch
      .map{[it[0]["library_id"]] + it}
      // add in singler results
      .join(singler_output_ch, by: 0, failOnMismatch: true, failOnDuplicate: true)
      // add in cell assign results
      .join(cellassign_output_ch, by: 0, failOnMismatch: true, failOnDuplicate: true)
      .map{it.drop(1)} // remove library_id
      .branch{
        // pull out libraries that actually have at least 1 type of annotation
        add_celltypes: (it[2].baseName != "NO_FILE") || (it[3].baseName != "NO_FILE")
        no_celltypes: true
      }


    // incorporate annotations into SCE object
    add_celltypes_to_sce(assignment_input_ch.add_celltypes)

    // mix in libraries without new celltypes
    // result is [meta, proccessed rds]
    celltyped_ch = assignment_input_ch.no_celltypes
      .map{[it[0], it[1]]}
      .mix(add_celltypes_to_sce.out)

    // add back in the unchanged sce files to the results
    export_channel = celltyped_ch
      .map{[it[0]["library_id"]] + it}
      // add in unfiltered and filtered sce files
      .join(sce_files_channel.map{[it[0]["library_id"], it[1], it[2]]},
            by: 0, failOnMismatch: true, failOnDuplicate: true)
      // rearrange to be [meta, unfiltered, filtered, processed]
      .map{library_id, meta, processed_sce, unfiltered_sce, filtered_sce ->
        [meta, unfiltered_sce, filtered_sce, processed_sce]
      }

  emit: export_channel

}
