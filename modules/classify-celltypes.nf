
process classify_singler {
    container params.SCPCATOOLS_CONTAINER
    publishDir (
      path: "${meta.celltype_publish_dir}",
      mode: 'copy',
      pattern: "${singler_dir}"
    )
    label 'mem_8'
    label 'cpus_4'
    input:
      tuple val(meta), path(processed_rds),
            path(singler_model_file),
            path(cellassign_reference_file), val(cellassign_reference_name)
    output:
      tuple val(meta), path(processed_rds),
            path(cellassign_reference_file), val(cellassign_reference_name),
            path(singler_dir)
    script:
      singler_dir = file(meta.singler_dir).name
      """
      # create output directory
      mkdir "${singler_dir}"

      classify_SingleR.R \
        --sce_file "${processed_rds}" \
        --singler_model_file "${singler_model_file}" \
        --output_singler_annotations_file "${singler_dir}/singler_annotations.tsv" \
        --output_singler_results_file "${singler_dir}/singler_results.rds" \
        --seed ${params.seed} \
        --threads ${task.cpus}

      # write out meta file
      echo "${Utils.makeJson(meta)}" > "${singler_dir}/scpca-meta.json"
      """
    stub:
      singler_dir = file(meta.singler_dir).name
      """
      mkdir "${singler_dir}"
      echo "${Utils.makeJson(meta)}" > "${singler_dir}/scpca-meta.json"
      """
}


process predict_cellassign {
  container params.SCPCATOOLS_CONTAINER
  publishDir "${params.checkpoints_dir}/celltype/${meta.library_id}", mode: 'copy'
  label 'mem_32'
  label 'cpus_12'
  input:
    tuple val(meta), path(processed_hdf5), path(cellassign_reference_mtx), val(ref_name)
  output:
    tuple val(meta), path(cellassign_predictions), val(ref_name)
  script:
    cellassign_predictions = "${meta.library_id}_predictions.tsv"
    """
    predict_cellassign.py \
      --input_hdf5_file ${processed_hdf5} \
      --output_predictions ${cellassign_predictions} \
      --reference ${cellassign_reference_mtx} \
      --seed ${params.seed} \
      --threads ${task.cpus}
    """
  stub:
    cellassign_predictions = "${meta.library_id}_predictions.tsv"
    """
    touch "${cellassign_predictions}"
    """
}

process classify_cellassign {
  container params.SCPCATOOLS_CONTAINER
  publishDir "${params.results_dir}/${meta.project_id}/${meta.sample_id}", mode: 'copy'
  label 'mem_4'
  label 'cpus_2'
  input:
    tuple val(meta), path(input_rds), path(cellassign_predictions), val(ref_name)
  output:
    tuple val(meta), path(annotated_rds)
  script:
    annotated_rds = "${meta.library_id}_annotated.rds"
    """
    classify_cellassign.R \
      --input_sce_file ${input_rds} \
      --output_sce_file ${annotated_rds} \
      --cellassign_predictions ${cellassign_predictions} \
      --reference_name ${ref_name}
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
         file(Utils.parseNA(it.singler_ref_file) ? "${params.singler_models_dir}/${it.singler_ref_file}" : empty_file),
         // cellassign reference file
         file(Utils.parseNA(it.cellassign_ref_file) ? "${params.cellassign_ref_dir}/${it.cellassign_ref_file}" : empty_file),
         // add ref name for cellassign since we cannot store it in the cellassign output
         // singler ref name does not need to be added because it is stored in the singler model
         Utils.parseNA(it.cellassign_ref_name)
        ]}

      // create input for typing: [meta, processed_sce, singler_model_file, cellassign_reference_file, cellassign_reference_name]
      celltype_input_ch = processed_sce_channel
        .map{[it[0]["project_id"]] + it}
        .combine(celltype_ch, by: 0)
        .map{it.drop(1)} // remove extra project ID
        // add celltype_publish_dir to the meta object for later use
        .map{
          it[0].celltype_publish_dir = "${params.checkpoints_dir}/celltype/${it[0].library_id}";
          it[0].singler_dir = "${it[0].celltype_publish_dir}/${it[0].library_id}_singler";
          it[0].cellassign_dir = "${it[0].celltype_publish_dir}/${it[0].library_id}_cellassign";
          it
        }

      // skip if no singleR model file
      singler_input_ch = celltype_input_ch
        .branch{
          missing_ref: it[2].name == "NO_FILE"
          do_singler: true
        }

      // perform singleR celltyping and export results
      classify_singler(singler_input_ch.do_singler)

      // singleR output channel
      singler_output_ch = singler_input_ch.missing_ref
        // drop singler model and add empty file for missing output (match classify_singler output)
        .map{meta, processed_sce, singler_model_file, cellassign_reference_file, cellassign_reference_name ->
             [meta, processed_sce, cellassign_reference_file, cellassign_reference_name, file(empty_file)]}
        // add in channel outputs
        .mix(classify_singler.out)


      // add back in the unchanged sce files
      // TODO update below with output channel results:
      export_channel = processed_sce_channel
        .map{[it[0]["library_id"]] + it}
        // add in unfiltered and filtered sce files
        .join(sce_files_channel.map{[it[0]["library_id"], it[1], it[2]]},
              by: 0, failOnMismatch: true, failOnDuplicate: true)
        // rearrange to be [meta, unfiltered, filtered, processed]
        .map{library_id, meta, processed_sce, unfiltered_sce, filtered_sce ->
            [meta, unfiltered_sce, filtered_sce, processed_sce]}

    emit: export_channel

}
