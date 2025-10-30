
include { getVersions; makeJson; readMeta; getMetaVal; parseNA } from '../lib/utils.nf'

process classify_singler {
  container Utils.pullthroughContainer(params.scpcatools_container, params.pullthrough_registry)
  publishDir (
    path: "${meta.celltype_checkpoints_dir}",
    mode: 'copy'
  )
  label 'mem_8'
  label 'cpus_4'
  tag "${meta.unique_id}"
  input:
    tuple val(meta), path(processed_rds), path(singler_model_file)
  output:
    tuple val(meta.unique_id), path(singler_dir)
  script:
    singler_dir = file(meta.singler_dir).name
    meta += getVersions(workflow, nextflow)
    meta_json = makeJson(meta)
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
    meta += getVersions(workflow, nextflow)
    meta_json = makeJson(meta)
    """
    mkdir "${singler_dir}"
    echo '${meta_json}' > "${singler_dir}/scpca-meta.json"
    touch "${singler_dir}/singler_results.rds"
    """
}


process classify_cellassign {
  container Utils.pullthroughContainer(params.scpcatools_scvi_container, params.pullthrough_registry)
    publishDir (
      path: "${meta.celltype_checkpoints_dir}",
      mode: 'copy'
    )
  label 'mem_max'
  label 'cpus_12'
  label 'long_running'
  tag "${meta.unique_id}"
  input:
    tuple val(meta), path(processed_rds), path(cellassign_reference_file)
  output:
    tuple val(meta.unique_id), path(cellassign_dir)
  script:
    cellassign_dir = file(meta.cellassign_dir).name
    meta += getVersions(workflow, nextflow)
    meta_json = makeJson(meta)
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
    meta += getVersions(workflow, nextflow)
    meta_json = makeJson(meta)
    """
    mkdir "${cellassign_dir}"
    echo '${meta_json}' > "${cellassign_dir}/scpca-meta.json"
    touch "${cellassign_dir}/cellassign_predictions.tsv"
    """
}

process classify_scimilarity {
  container Utils.pullthroughContainer(params.scpcatools_scimilarity_container, params.pullthrough_registry)
    publishDir (
      path: "${meta.celltype_checkpoints_dir}",
      mode: 'copy'
    )
  label 'mem_96'
  label 'cpus_4'
  tag "${meta.unique_id}"
  input:
    tuple val(meta), path(processed_rds), path(scimilarity_model_dir), path(scimilarity_ontology_map_file)
  output:
    tuple val(meta.unique_id), path(scimilarity_dir)
  script:
    scimilarity_dir = file(meta.scimilarity_dir).name
    def meta_json = makeJson(meta + getVersions(workflow, nextflow))
    """
    # create output directory
    mkdir "${scimilarity_dir}"

    # Convert SCE to AnnData
    sce_to_anndata.R \
        --input_sce_file "${processed_rds}" \
        --output_rna_h5 "processed.h5ad"

    # only run cell assign if a h5ad file was able to be created
    # otherwise create an empty predictions file
    if [[ -e "processed.h5ad" ]]; then

      # Run SCimilarity
      run_scimilarity.py \
        --model_dir "${scimilarity_model_dir}" \
        --processed_h5ad_file "processed.h5ad" \
        --ontology_map_file ${scimilarity_ontology_map_file} \
        --predictions_tsv "${scimilarity_dir}/scimilarity_predictions.tsv" \
        --seed ${params.seed}

    else
      touch "${scimilarity_dir}/scimilarity_predictions.tsv"
    fi

    # write out meta file
    echo '${meta_json}' > "${scimilarity_dir}/scpca-meta.json"

    """
  stub:
    scimilarity_dir = file(meta.scimilarity_dir).name
    def meta_json = makeJson(meta + getVersions(workflow, nextflow))
    """
    mkdir "${scimilarity_dir}"
    echo '${meta_json}' > "${scimilarity_dir}/scpca-meta.json"
    touch "${scimilarity_dir}/scimilarity_predictions.tsv"
    """
}

process add_celltypes_to_sce {
  container Utils.pullthroughContainer(params.scpcatools_slim_container, params.pullthrough_registry)
  label 'mem_4'
  label 'cpus_2'
  tag "${meta.unique_id}"
  input:
    tuple val(meta), path(processed_rds), path(singler_dir), path(cellassign_dir), path(scimilarity_dir)
    path(celltype_ref_metadata) // TSV file of references metadata needed for CellAssign only
    path(panglao_ref_file) // used for assigning ontology IDs for CellAssign results
    path(consensus_ref_file) // used for assigning consensus cell types if both SingleR and CellAssign are used
    path(validation_ref_file) // maps consensus cell types to cell type groups, for counting normal reference cells
    path(diagnosis_celltypes_file) // maps broad diagnoses to cell type groups, for counting normal reference cells
    path(diagnosis_groups_file) // maps broad diagnoses to cell type groups, for counting normal reference cells
  output:
    tuple val(meta), path(annotated_rds), env("REFERENCE_CELL_COUNT"), env("REFERENCE_CELL_HASH")
  script:
    annotated_rds = "${meta.unique_id}_processed_annotated.rds"
    def singler_results = singler_dir ? "${singler_dir}/singler_results.rds": ""
    def cellassign_predictions = cellassign_dir ? "${cellassign_dir}/cellassign_predictions.tsv" : ""
    def scimilarity_results = scimilarity_dir ? "${scimilarity_dir}/scimilarity_predictions.tsv" : ""

    """
    add_celltypes_to_sce.R \
      --input_sce_file ${processed_rds} \
      --output_sce_file ${annotated_rds} \
      --singler_results  "${singler_results}" \
      --singler_model_file "${meta.singler_model_file}" \
      --cellassign_predictions  "${cellassign_predictions}" \
      --cellassign_ref_file "${meta.cellassign_reference_file}" \
      --celltype_ref_metafile "${celltype_ref_metadata}" \
      --panglao_ontology_ref "${panglao_ref_file}" \
      --scimilarity_results "${scimilarity_results}" \
      --scimilarity_model_dir "${meta.scimilarity_model_dir}" \
      --consensus_celltype_ref "${consensus_ref_file}" \
      --consensus_validation_ref "${validation_ref_file}" \
      --diagnosis_celltype_ref "${diagnosis_celltypes_file}" \
      --diagnosis_groups_ref "${diagnosis_groups_file}" \
      --reference_cell_count_file "reference_cell_count.txt" \
      --reference_cell_hash_file "reference_cell_hash.txt" \

      # save so we can export as environment variable
      REFERENCE_CELL_COUNT=\$(cat "reference_cell_count.txt")
      REFERENCE_CELL_HASH=\$(cat "reference_cell_hash.txt")
    """
  stub:
    annotated_rds = "${meta.unique_id}_processed_annotated.rds"
    """
    touch ${annotated_rds}

    # Set to a value guaranteed to pass the threshold and run
    REFERENCE_CELL_COUNT=${params.infercnv_min_reference_cells + 1}
    REFERENCE_CELL_HASH=""
    """
}


workflow annotate_celltypes {
  take: sce_files_channel // channel of meta, unfiltered_sce, filtered_sce, processed_sce
  main:

    // read in sample metadata and make a list of cell line samples; these won't be cell typed
    cell_line_samples = channel.fromPath(params.sample_metafile)
      .splitCsv(header: true, sep: '\t')
      .filter{ it.is_cell_line.toBoolean() }
      .map{ it -> it.scpca_sample_id }
      .toList()

    // branch to cell type the non-cell line libraries only
    sce_files_channel_branched = sce_files_channel
      .branch{ meta, _unfiltered, _filtered, _processed ->
        cell_line: meta.sample_id.split(",").every{ it in cell_line_samples.getVal() }
        // only run cell typing on tissue samples
        tissue: true
      }

    // get just the meta and processed sce from the tissue (not cell line) samples
    processed_sce_channel = sce_files_channel_branched.tissue
      .map{ meta, _unfiltered, _filtered, processed ->
        [meta, processed]
      }

    // channel with [project_id, singler_model_file, cellassign_reference_file]
    celltype_ch = channel.fromPath(params.project_celltype_metafile)
      .splitCsv(header: true, sep: '\t')
      .map{ it ->
        def singler_model_file = parseNA(it.singler_ref_file) ? "${params.singler_models_dir}/${it.singler_ref_file}" : ''
        def cellassign_ref_file = parseNA(it.cellassign_ref_file) ? "${params.cellassign_ref_dir}/${it.cellassign_ref_file}" : ''

        [it.scpca_project_id, singler_model_file, cellassign_ref_file]
      }

    // create input for typing: [augmented meta, processed_sce]
    celltype_input_ch = processed_sce_channel
      .map{ meta, processed_sce ->
        [meta.project_id, meta, processed_sce]
      }
      .combine(celltype_ch, by: 0)
      // current contents: [project_id, meta, processed_sce, singler_model_file, cellassign_reference_file]
      // add values to meta for later use
      .map{ _project_id, meta_in, processed_sce, singler_model_file, cellassign_reference_file ->
        def meta = meta_in.clone() // local copy for safe modification
        // results directories
        meta.celltype_checkpoints_dir = "${params.checkpoints_dir}/celltype/${meta.library_id}"
        meta.singler_dir = "${meta.celltype_checkpoints_dir}/${meta.unique_id}_singler"
        meta.cellassign_dir = "${meta.celltype_checkpoints_dir}/${meta.unique_id}_cellassign"
        meta.scimilarity_dir = "${meta.celltype_checkpoints_dir}/${meta.unique_id}_scimilarity"
        // reference files
        meta.singler_model_file = singler_model_file
        meta.cellassign_reference_file = cellassign_reference_file
        meta.scimilarity_model_dir = params.scimilarity_model_dir
        meta.scimilarity_ontology_map_file = params.scimilarity_ontology_map_file
        // output files
        meta.singler_results_file = "${meta.singler_dir}/singler_results.rds"
        meta.cellassign_predictions_file = "${meta.cellassign_dir}/cellassign_predictions.tsv"
        meta.scimilarity_predictions_file = "${meta.scimilarity_dir}/scimilarity_predictions.tsv"

        // return simplified input:
        [meta, processed_sce]
      }

    /////////////////////////////////////////////////////
    //                  SingleR                        //
    /////////////////////////////////////////////////////

    // creates [meta, processed sce, singler model file]
    singler_input_ch = celltype_input_ch
      // add in singler model file
      .map{ meta, processed_sce ->
        def singler_model = meta.singler_model_file ? file(meta.singler_model_file, checkIfExists: true) : []
        [meta, processed_sce, singler_model]
      }
      // skip if no singleR model file or if singleR results are already present
      .branch{ meta, _processed_sce, singler_model ->
        def stored_singler_model_file = getMetaVal(file("${meta.singler_dir}/scpca-meta.json"), "singler_model_file")
        skip_singler: (
          !params.repeat_celltyping
          && file(meta.singler_results_file).exists()
          && meta.singler_model_file == stored_singler_model_file
        )
        missing_ref: singler_model == []
        do_singler: true
      }


    // perform singleR celltyping and export results
    classify_singler(singler_input_ch.do_singler)

    // singleR output channel: [unique id, singler_results]
    singler_output_ch = singler_input_ch.skip_singler
      // provide existing singler results dir for those we skipped
      .map{ meta, _processed_sce, _singler_model ->
        [
          meta.unique_id,
          file(meta.singler_dir, type: 'dir', checkIfExists: true)
        ]
      }
      // add in missing ref samples
      .mix(singler_input_ch.missing_ref.map{ it -> [it[0].unique_id, []] })
      // add in channel outputs
      .mix(classify_singler.out)


    /////////////////////////////////////////////////////
    //                  CellAssign                     //
    /////////////////////////////////////////////////////

    // create cellassign input channel: [meta, processed sce, cellassign reference file]
    cellassign_input_ch = celltype_input_ch
      // add in cellassign reference
      .map{ meta, processed_sce ->
        def cellassign_ref = meta.cellassign_reference_file ? file(meta.cellassign_reference_file, checkIfExists: true) : []
        [meta, processed_sce, cellassign_ref]
      }
      // skip if no cellassign reference file or reference name is not defined
      .branch{ meta, _processed_sce, cellassign_ref ->
        def stored_cellassign_reference_file = getMetaVal(file("${meta.cellassign_dir}/scpca-meta.json"), "cellassign_reference_file")
        skip_cellassign: (
          !params.repeat_celltyping
          && file(meta.cellassign_predictions_file).exists()
          && meta.cellassign_reference_file == stored_cellassign_reference_file
        )
        missing_ref: cellassign_ref == []
        do_cellassign: true
      }


    // perform CellAssign celltyping and export results
    classify_cellassign(cellassign_input_ch.do_cellassign)

    // cellassign output channel: [unique id, cellassign_dir]
    cellassign_output_ch = cellassign_input_ch.skip_cellassign
      // provide existing cellassign predictions dir for those we skipped
      .map{ meta, _processed_sce, _cellassign_ref ->
        [
          meta.unique_id,
          file(meta.cellassign_dir, type: 'dir', checkIfExists: true)
        ]
      }
      // add missing ref samples
      .mix(cellassign_input_ch.missing_ref.map{ it -> [it[0].unique_id, []] })
      // add in channel outputs
      .mix(classify_cellassign.out)


    /////////////////////////////////////////////////////
    //                  SCimilarity                    //
    /////////////////////////////////////////////////////

    // create scimilarity input channel: [meta, processed sce, scimilarity dir, ontology map]
    scimilarity_input_ch = celltype_input_ch
      // add in cellassign references
      .map{ meta, processed_sce ->
        [
          meta,
          processed_sce,
          file(meta.scimilarity_model_dir, type: 'dir', checkIfExists: true),
          file(meta.scimilarity_ontology_map_file, checkIfExists: true)
        ]
      }
      // skip if scimilarity results exist and path to the model has not changed
      .branch{ it ->
        def stored_scimilarity_model_dir = getMetaVal(file("${it[0].scimilarity_dir}/scpca-meta.json"), "scimilarity_model_dir")

        skip_scimilarity: (
          !params.repeat_celltyping
          && file(it[0].scimilarity_predictions_file).exists()
          && it[0].scimilarity_model_dir == stored_scimilarity_model_dir
        )
        do_scimilarity: true
      }


    // run SCimilarity and export results
    classify_scimilarity(scimilarity_input_ch.do_scimilarity)

    // scimilarity output channel: [unique id, scimilarity_dir]
    scimilarity_output_ch = scimilarity_input_ch.skip_scimilarity
      // provide existing scimilarity predictions dir for those we skipped
      .map{ meta, _processed_sce, _scimilarity_model_dir, _scimilarity_ontology_map_file ->
        [
          meta.unique_id,
          file(meta.scimilarity_dir, type: 'dir', checkIfExists: true)
        ]
      }
      // add in channel outputs
      .mix(classify_scimilarity.out)

    /////////////////////////////////////////////////////
    //               Add celltypes to object           //
    /////////////////////////////////////////////////////

    // prepare input for process to add celltypes to the processed SCE
    // result is [meta, processed rds, singler dir, cellassign dir, scimilarity dir]
    assignment_input_ch = celltype_input_ch
      .map{ meta, processed_sce ->
        [meta.unique_id, meta, processed_sce]
      }
      // add in singler results
      .join(singler_output_ch, by: 0, failOnMismatch: true, failOnDuplicate: true)
      // add in cell assign results
      .join(cellassign_output_ch, by: 0, failOnMismatch: true, failOnDuplicate: true)
      // add in scimilarity results
      .join(scimilarity_output_ch, by: 0, failOnMismatch: true, failOnDuplicate: true)
      .map{ it -> it.drop(1) } // remove unique id
      // pull out libraries that actually have at least 1 type of annotation

      .branch{ _meta, _processed, singler_dir, cellassign_dir, scimilarity_dir ->
        def has_annotation = [singler_dir, cellassign_dir, scimilarity_dir].any()
        add_celltypes: has_annotation
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
      params.diagnosis_celltypes_file ? file(params.diagnosis_celltypes_file, checkIfExists: true) : [], // maps broad diagnoses to cell type groups, for counting normal reference cells
      params.diagnosis_groups_file ? file(params.diagnosis_groups_file, checkIfExists: true) : [] // maps sample diagnoses to broad diagnoses, for counting normal reference cells
    )

    // add inferCNV logic to meta
    added_celltypes_ch = add_celltypes_to_sce.out
      .map{ meta_in, annotated_sce, cell_count, cell_hash ->
        def meta = meta_in.clone() // local copy for safe modification
        // ensure the count is saved as an integer: either the integer value, or null if it was an
        // empty string since we can do future math comparisons with null
        meta.infercnv_reference_cell_count = cell_count ? cell_count.toInteger() : null
        meta.infercnv_reference_cell_hash = cell_hash

        // return only meta and annotated_sce
        [meta, annotated_sce]
      }


    // mix in libraries without new celltypes
    // result is [meta, processed rds]
    celltyped_ch = assignment_input_ch.no_celltypes
      .map{ meta, processed, _singler_dir, _cellassign_dir, _scimilarity_dir ->
        [meta, processed]
      }
      .mix(added_celltypes_ch)


    // add back in the unchanged sce files to the results
    export_channel = celltyped_ch
      .map{ meta, processed_sce ->
        [meta.unique_id, meta, processed_sce]
      }
      // add in unfiltered and filtered sce files, for tissue samples only
      .join(
        sce_files_channel_branched.tissue
          .map{ meta, unfiltered, filtered, _processed ->
            [meta.unique_id, unfiltered, filtered]
          },
        by: 0, failOnMismatch: true, failOnDuplicate: true
      )
      // rearrange to be [meta, unfiltered, filtered, processed]
      .map{ _unique_id, meta, processed_sce, unfiltered_sce, filtered_sce ->
        [meta, unfiltered_sce, filtered_sce, processed_sce]
      }
      // mix in cell line libraries which were not cell typed
      .mix(sce_files_channel_branched.cell_line)

  emit: export_channel

}
