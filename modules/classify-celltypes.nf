
process classify_singleR {
    container params.SCPCATOOLS_CONTAINER
    label 'mem_8'
    label 'cpus_4'
    publishDir "${params.results_dir}/${meta.project_id}/${meta.sample_id}", mode: 'copy'
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
  label 'mem_32'
  label 'cpus_12'
  input:
    tuple val(meta), path(processed_hdf5), path(cellassign_reference_mtx)
  output:
    tuple val(meta), path(cellassign_predictions)
  script:
    cellassign_predictions = "${meta.library_id}_predictions.tsv"
    """
    classify_cellassign.py \
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

workflow annotate_celltypes {
    take: processed_sce_channel
    main:
      // channel with celltype model and project ids
      celltype_ch = Channel.fromPath(params.celltype_project_metafile)
        .splitCsv(header: true, sep: '\t')
        .map{[
          project_id = it.scpca_project_id,
          singler_model_file = "${params.singler_models_dir}/${it.singler_ref_file}",
          cellassign_ref_file = "${params.cellassign_ref_dir}/${it.cellassign_ref_file}"
        ]}

      // create channel grouped_celltype_ch as: [meta, processed sce object, SingleR reference model]
      // input processed_sce_channel is [meta, unfiltered, filtered, processed]
      grouped_celltype_ch = processed_sce_channel
        .map{[it[0]["project_id"]] + it}
        .combine(celltype_ch, by: 0)
        .map{it.drop(1)} // remove extra project ID

      // creates [meta, processed, SingleR reference model]
      singler_input_ch = grouped_celltype_ch
        .map{meta, processed_rds, processed_hdf5, singler_model_file, cellassign_ref_file -> tuple(meta,
                                                                                                  processed_rds,
                                                                                                  singler_model_file
                                                                                                  )}

      cellassign_input_ch = grouped_celltype_ch
        .map{meta, processed_rds, processed_hdf5, singler_model_file, cellassign_ref_file -> tuple(meta,
                                                                                                  processed_hdf5,
                                                                                                  cellassign_ref_file
                                                                                                  )}


      classify_singleR(singler_input_ch)

      predict_cellassign(cellassign_input_ch)

    emit: tuple(classify_singleR.out, predict_cellassign.out)

}
