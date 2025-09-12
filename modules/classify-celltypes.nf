
process classify_singler {
  container params.SCPCATOOLS_CONTAINER
  publishDir (
    path: "${meta.celltype_checkpoints_dir}",
    mode: 'copy',
    pattern: "${singler_dir}"
  )
  label 'mem_8'
  label 'cpus_4'
  tag "${meta.library_id}"
  input:
    tuple val(meta), path(processed_rds), path(singler_model_file)
  output:
    tuple val(meta.unique_id), path(singler_dir)
  script:
    singler_dir = file(meta.singler_dir).name
    meta += Utils.getVersions(workflow, nextflow)
    meta_json = Utils.makeJson(meta)
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
    echo '${meta_json}' > "${singler_dir}/scpca-meta.json"
    """
  stub:
    singler_dir = file(meta.singler_dir).name
    meta += Utils.getVersions(workflow, nextflow)
    meta_json = Utils.makeJson(meta)
    """
    mkdir "${singler_dir}"
    echo '${meta_json}' > "${singler_dir}/scpca-meta.json"
    touch "${singler_dir}/singler_results.rds"
    """
}


process classify_cellassign {
  container params.SCPCATOOLS_SCVI_CONTAINER
    publishDir (
      path: "${meta.celltype_checkpoints_dir}",
      mode: 'copy',
      pattern: "${cellassign_dir}"
    )
  label 'mem_max'
  label 'cpus_12'
  label 'long_running'
  tag "${meta.library_id}"
  input:
    tuple val(meta), path(processed_rds), path(cellassign_reference_file)
  output:
    tuple val(meta.unique_id), path(cellassign_dir)
  script:
    cellassign_dir = file(meta.cellassign_dir).name
    meta += Utils.getVersions(workflow, nextflow)
    meta_json = Utils.makeJson(meta)
    """
    # create output directory
    mkdir "${cellassign_dir}"

    # Convert SCE to AnnData
    sce_to_anndata.R \
        --input_sce_file "${processed_rds}" \
        --output_rna_h5 "processed.h5ad"

    # only run cell assign if a h5ad file was able to be created
    # otherwise create an empty predictions file
    if [[ -e "processed.h5ad" ]]; then

      # Run CellAssign
      predict_cellassign.py \
        --anndata_file "processed.h5ad" \
        --output_predictions "${cellassign_dir}/cellassign_predictions.tsv" \
        --reference "${cellassign_reference_file}" \
        --seed ${params.seed} \
        --threads ${task.cpus}

    else
      touch "${cellassign_dir}/cellassign_predictions.tsv"
    fi

    # write out meta file
      echo '${meta_json}' > "${cellassign_dir}/scpca-meta.json"

    """
  stub:
    cellassign_dir = file(meta.cellassign_dir).name
    meta += Utils.getVersions(workflow, nextflow)
    meta_json = Utils.makeJson(meta)
    """
    mkdir "${cellassign_dir}"
    echo '${meta_json}' > "${cellassign_dir}/scpca-meta.json"
    touch "${cellassign_dir}/cellassign_predictions.tsv"
    """
}

process add_celltypes_to_sce {
  container params.SCPCATOOLS_SLIM_CONTAINER
  label 'mem_4'
  label 'cpus_2'
  tag "${meta.library_id}"
  input:
    tuple val(meta), path(processed_rds), path(singler_dir), path(cellassign_dir)
    path(celltype_ref_metadata) // TSV file of references metadata needed for CellAssign only
    path(panglao_ref_file) // used for assigning ontology IDs for CellAssign results
    path(consensus_ref_file) // used for assigning consensus cell types if both SingleR and CellAssign are used
    path(validation_ref_file) // maps consensus cell types to cell type groups, for counting normal reference cells
    path(diagnosis_celltypes_file, name: 'diagnosis_celltypes.txt') // maps broad diagnoses to cell type groups, for counting normal reference cells
    path(diagnosis_groups_file, name: 'diagnosis_groups.txt') // maps broad diagnoses to cell type groups, for counting normal reference cells
  output:
    tuple val(meta), path(annotated_rds), env("REFERENCE_CELL_COUNT")
  script:
    annotated_rds = "${meta.library_id}_processed_annotated.rds"
    singler_present = "${singler_dir.name}" != "NO_FILE"
    singler_results = "${singler_dir}/singler_results.rds"
    cellassign_present = "${cellassign_dir.name}" != "NO_FILE"
    cellassign_predictions = "${cellassign_dir}/cellassign_predictions.tsv"
    """
    add_celltypes_to_sce.R \
      --input_sce_file ${processed_rds} \
      --output_sce_file ${annotated_rds} \
      ${singler_present ? "--singler_results  ${singler_results}" : ''} \
      ${singler_present ? "--singler_model_file ${meta.singler_model_file}" : ''} \
      ${cellassign_present ? "--cellassign_predictions  ${cellassign_predictions}" : ''} \
      ${cellassign_present ? "--cellassign_ref_file ${meta.cellassign_reference_file}" : ''} \
      ${cellassign_present ? "--celltype_ref_metafile ${celltype_ref_metadata}" : ''} \
      ${cellassign_present ? "--panglao_ontology_ref ${panglao_ref_file}" : ''} \
      --consensus_celltype_ref ${consensus_ref_file} \
      --consensus_validation_ref ${validation_ref_file} \
      --diagnosis_celltype_ref "diagnosis_celltypes.txt" \
      --diagnosis_groups_ref "diagnosis_groups.txt" \
      --reference_cell_count_file "reference_cell_count.txt"

      # save so we can export as environment variable
      REFERENCE_CELL_COUNT=\$(cat "reference_cell_count.txt")
    """
  stub:
    annotated_rds = "${meta.library_id}_processed_annotated.rds"
    """
    touch ${annotated_rds}

    # Set to a value guaranteed to pass the threshold
    REFERENCE_CELL_COUNT=${params.infercnv_min_reference_cells + 1}
    """
}


workflow annotate_celltypes {
  take: sce_files_channel // channel of meta, unfiltered_sce, filtered_sce, processed_sce
  main:
    def empty_file = "${projectDir}/assets/NO_FILE"
    // read in sample metadata and make a list of cell line samples; these won't be cell typed
    cell_line_samples = Channel.fromPath(params.sample_metafile)
      .splitCsv(header: true, sep: '\t')
      .filter{it.is_cell_line.toBoolean()}
      .map{it.scpca_sample_id}
      .toList()

    // branch to cell type the non-cell line libraries only
    sce_files_channel_branched = sce_files_channel
     .branch{
        cell_line: it[0]["sample_id"].split(",").every{it in cell_line_samples.getVal()}
        // only run cell typing on tissue samples
        tissue: true
      }

    // get just the meta and processed sce from the tissue (not cell line) samples
    processed_sce_channel = sce_files_channel_branched.tissue.map{[it[0], it[3]]}

    // channel with celltype model and project ids
    celltype_ch = Channel.fromPath(params.project_celltype_metafile)
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
      .map{ meta, processed_sce -> tuple(
        meta.project_id,
        meta,
        processed_sce
        )}
      .combine(celltype_ch, by: 0)
      // current contents: [project_id, meta, processed_sce, singler_model_file, cellassign_reference_file]
      // add values to meta for later use
      .map{ _project_id, meta_in, processed_sce, singler_model_file, cellassign_reference_file ->
        def meta = meta_in.clone(); // local copy for safe modification
        meta.celltype_checkpoints_dir = "${params.checkpoints_dir}/celltype/${meta.library_id}";
        meta.singler_dir = "${meta.celltype_checkpoints_dir}/${meta.unique_id}_singler";
        meta.cellassign_dir = "${meta.celltype_checkpoints_dir}/${meta.unique_id}_cellassign";
        meta.singler_model_file = singler_model_file;
        meta.cellassign_reference_file = cellassign_reference_file;
        meta.singler_results_file = "${meta.singler_dir}/singler_results.rds";
        meta.cellassign_predictions_file = "${meta.cellassign_dir}/cellassign_predictions.tsv";
        // return simplified input:
        [meta, processed_sce]
      }


    // creates [meta, processed sce, singler model file]
    singler_input_ch = celltype_input_ch
      // add in singler model or empty file
      .map{it.toList() + [file(it[0].singler_model_file ?: empty_file, checkIfExists: true)]}
      // skip if no singleR model file or if singleR results are already present
      .branch{
        skip_singler: (
          !params.repeat_celltyping
          && file(it[0].singler_results_file).exists()
          && Utils.getMetaVal(file("${it[0].singler_dir}/scpca-meta.json"), "singler_model_file") == "${it[0].singler_model_file}"
        )
        missing_ref: it[2].name == "NO_FILE"
        do_singler: true
      }


    // perform singleR celltyping and export results
    classify_singler(singler_input_ch.do_singler)

    // singleR output channel: [unique id, singler_results]
    singler_output_ch = singler_input_ch.skip_singler
      // provide existing singler results dir for those we skipped
      .map{ meta, processed_sce, singler_model -> tuple(
        meta.unique_id,
        file(meta.singler_dir, type: 'dir', checkIfExists: true)
      )}
      // add empty file for missing ref samples
      .mix(singler_input_ch.missing_ref.map{[it[0]["unique_id"], file(empty_file, checkIfExists: true)]} )
      // add in channel outputs
      .mix(classify_singler.out)

    // create cellassign input channel: [meta, processed sce, cellassign reference file]
    cellassign_input_ch = celltype_input_ch
      // add in cellassign reference
      .map{it.toList() + [file(it[0].cellassign_reference_file ?: empty_file, checkIfExists: true)]}
      // skip if no cellassign reference file or reference name is not defined
      .branch{
        skip_cellassign: (
          !params.repeat_celltyping
          && file(it[0].cellassign_predictions_file).exists()
          && Utils.getMetaVal(file("${it[0].cellassign_dir}/scpca-meta.json"), "cellassign_reference_file") == "${it[0].cellassign_reference_file}"
        )
        missing_ref: it[2].name == "NO_FILE"
        do_cellassign: true
      }


    // perform CellAssign celltyping and export results
    classify_cellassign(cellassign_input_ch.do_cellassign)

    // cellassign output channel: [unique id, cellassign_dir]
    cellassign_output_ch = cellassign_input_ch.skip_cellassign
      // provide existing cellassign predictions dir for those we skipped
      .map{ meta, processed_sce, cellassign_ref -> tuple(
        meta.unique_id,
        file(meta.cellassign_dir, type: 'dir', checkIfExists: true)
      )}
      // add empty file for missing ref samples
      .mix(cellassign_input_ch.missing_ref.map{[it[0]["unique_id"], file(empty_file, checkIfExists: true)]} )
      // add in channel outputs
      .mix(classify_cellassign.out)

    // prepare input for process to add celltypes to the processed SCE
    // result is [meta, processed rds, singler dir, cellassign dir]
    assignment_input_ch = celltype_input_ch
      .map{ meta, processed_sce -> tuple(
        meta.unique_id,
        meta,
        processed_sce
      )}
      // add in singler results
      .join(singler_output_ch, by: 0, failOnMismatch: true, failOnDuplicate: true)
      // add in cell assign results
      .join(cellassign_output_ch, by: 0, failOnMismatch: true, failOnDuplicate: true)
      .map{it.drop(1)} // remove unique id
      .branch{
        // pull out libraries that actually have at least 1 type of annotation
        add_celltypes: (it[2].baseName != "NO_FILE") || (it[3].baseName != "NO_FILE")
        no_celltypes: true
      }


    // incorporate annotations into SCE object
    // outputs [meta, annotated processed rds, reference cell count]
    add_celltypes_to_sce(
      assignment_input_ch.add_celltypes,
      file(params.celltype_ref_metadata), // file with CellAssign reference organs
      file(params.panglao_ref_file), // used for assigning ontology IDs for CellAssign results
      file(params.consensus_ref_file), // used for assigning consensus cell types if both SingleR and CellAssign are used
      file(params.validation_groups_file),  // maps consensus cell types to cell type groups, for counting normal reference cells
      file(params.diagnosis_celltypes_file ?: empty_file, checkIfExists: true), // maps broad diagnoses to cell type groups, for counting normal reference cells
      file(params.diagnosis_groups_file ?: empty_file, checkIfExists: true) // maps sample diagnoses to broad diagnoses, for counting normal reference cells
    )

    // add inferCNV logic to meta
    added_celltypes_ch = add_celltypes_to_sce.out
      .map{ meta_in, annotated_sce, cell_count ->
        def meta = meta_in.clone(); // local copy for safe modification
        // ensure it's saved as an integer: either the integer value, or null if it was NA
        meta.infercnv_reference_cell_count = Utils.parseNA(cell_count) == "" ? null : cell_count.toInteger();
        // return only meta and annotated_sce
        [meta, annotated_sce]
      }


    // mix in libraries without new celltypes
    // result is [meta, processed rds]
    celltyped_ch = assignment_input_ch.no_celltypes
      .map{[it[0], it[1]]}
      .mix(added_celltypes_ch)

    // add back in the unchanged sce files to the results
    export_channel = celltyped_ch
      .map{meta, processed_sce -> tuple(
        meta.unique_id,
        meta,
        processed_sce
        )}
      // add in unfiltered and filtered sce files, for tissue samples only
      .join(
        sce_files_channel_branched.tissue.map{[it[0]["unique_id"], it[1], it[2]]},
        by: 0, failOnMismatch: true, failOnDuplicate: true
      )
      // rearrange to be [meta, unfiltered, filtered, processed]
      .map{_unique_id, meta, processed_sce, unfiltered_sce, filtered_sce ->
        [meta, unfiltered_sce, filtered_sce, processed_sce]
      }
      // mix in cell line libraries which were not cell typed
      .mix(sce_files_channel_branched.cell_line)

  emit: export_channel

}
