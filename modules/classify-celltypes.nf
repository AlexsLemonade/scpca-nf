
process classify_singleR {
    container params.SCPCATOOLS_CONTAINER
    label 'mem_8'
    label 'cpus_4'
    input:
        tuple val(meta), path(processed_rds), path(singler_model_file)
    output:
        tuple val(meta), path(annotated_rds)
    script:
      annotated_rds = "${meta.library_id}_annotated.rds"
      """
      classify_SingleR.R \
        --input_sce_file ${processed_rds} \
        --output_sce_file ${annotated_rds} \
        --singler_model_file ${singler_model_file} \
        --label_name ${params.label_name} \
        --seed ${params.seed} \
        --threads ${task.cpus}
      """
    stub:
      annotated_rds = "${meta.library_id}_annotated.rds"
      """
      touch "${annotated_rds}"
      """
}

process predict_cellassign {
  container params.SCPCATOOLS_CONTAINER
  publishDir "${params.results_dir}/${meta.project_id}/${meta.sample_id}", mode: 'copy'
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

workflow annotate_celltypes {
    take: processed_sce_channel
    main:
      // channel with celltype model and project ids
      celltype_ch = Channel.fromPath(params.celltype_project_metafile)
        .splitCsv(header: true, sep: '\t')
        .map{[
          project_id = it.scpca_project_id,
          singler_model_file = "${params.singler_models_dir}/${it.singler_ref_file}",
          cellassign_ref_file = "${params.cellassign_ref_dir}/${it.cellassign_ref_file}",
          cellassign_ref_name = it.cellassign_ref_name
        ]}

      // create channel grouped_celltype_ch as: [meta, processed sce object, SingleR reference model]
      // input processed_sce_channel is [meta, unfiltered, filtered, processed]
      grouped_celltype_ch = processed_sce_channel
        .map{[it[0]["project_id"]] + it}
        .combine(celltype_ch, by: 0)
        .map{it.drop(1)} // remove extra project ID

      // creates [meta, processed, SingleR reference model]
      singler_input_ch = grouped_celltype_ch
        .map{meta, processed_rds, processed_hdf5, singler_model_file, cellassign_ref_file, cellassign_ref_name -> tuple(meta,
                                                                                                                        processed_rds,
                                                                                                                        singler_model_file
                                                                                                                        )}

      // creates [meta, processed hdf5, cellassign ref file, cell assign ref name]
      cellassign_input_ch = grouped_celltype_ch
        .map{meta, processed_rds, processed_hdf5, singler_model_file, cellassign_ref_file, cellassign_ref_name -> tuple(meta,
                                                                                                                        processed_hdf5,
                                                                                                                        cellassign_ref_file,
                                                                                                                        cellassign_ref_name
                                                                                                                        )}

      // get SingleR cell type assignments and add them to SCE
      classify_singleR(singler_input_ch)

      // get cellassign predictions file
      predict_cellassign(cellassign_input_ch)

      // add cellassign annotations to the object with singleR results
      all_celltype_assignments_ch = classify_singleR.out
        // combines using meta from both singleR and cellassign, they should be the same
        // resulting tuple should be [meta, singleR annotated rds, cellassign predictions]
        .combine(predict_cellassign.out, by: 0)

      classify_cellassign(all_celltype_assignments_ch)

    emit: classify_cellassign.out

}
