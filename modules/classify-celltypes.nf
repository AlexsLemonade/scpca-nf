
process classify_singleR {
    container params.SCPCATOOLS_CONTAINER
    label 'mem_8'
    label 'cpus_4'
    input:
        tuple val(meta), path(unfiltered_rds), path(filtered_rds), path(processed_rds), path(singler_model_file)
    output:
        tuple val(meta), path(unfiltered_rds), path(filtered_rds), path(processed_rds)
    script:
      meta["has_singler"] = true
      """
      classify_SingleR.R \
        --sce_file ${processed_rds} \
        --singler_model_file ${singler_model_file} \
        --label_name ${params.singler_label_name} \
        --seed ${params.seed} \
        --threads ${task.cpus}
      """
  stub:
  """
  # nothing to do since files don't move
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
    tuple val(meta), path(unfiltered_rds), path(filtered_rds), path(processed_rds), path(cellassign_predictions), val(ref_name)
  output:
    tuple val(meta), path(unfiltered_rds), path(filtered_rds), path(processed_rds)
  script:
    meta["has_cellassign"] = true
    """
    classify_cellassign.R \
      --sce_file ${processed_rds} \
      --cellassign_predictions ${cellassign_predictions} \
      --reference_name ${ref_name}
    """
}

workflow annotate_celltypes {
    take: celltype_ch
    main:
      
      // channel with celltype model and project ids
      celltype_model_ch = Channel.fromPath(params.celltype_project_metafile)
        .splitCsv(header: true, sep: '\t')
        .map{[
          project_id = it.scpca_project_id,
          singler_model_file = "${params.singler_models_dir}/${it.singler_ref_file}",
          cellassign_ref_file = "${params.cellassign_ref_dir}/${it.cellassign_ref_file}",
          // add ref name for cellassign since we cannot store it in the cellassign output
          // note that singler ref name doesn't need to be added as it is stored in the singler model
          cellassign_ref_name = it.cellassign_ref_name
        ]}

      // From input, create channel with grouped meta, processed sce object, and all references to use
      // input celltype_ch is [[meta], processed_hdf5, unfiltered_rds, filtered_rds, processed_rds]
      grouped_celltype_ch = celltype_ch
        .map{[it[0]["project_id"]] + it}
        .combine(celltype_model_ch, by: 0)
        .map{it.drop(1)} // remove extra project ID
        .groupTuple(by: 0) // group by meta

      // creates input for singleR [meta, unfiltered, filtered, processed, SingleR reference model]
      singler_input_ch = grouped_celltype_ch
        .map{meta, processed_hdf5, unfiltered_rds, filtered_rds, processed_rds, singler_model_file, cellassign_ref_file, cellassign_ref_name -> tuple(meta,
                                                                                                                                                      unfiltered_rds,
                                                                                                                                                      filtered_rds,
                                                                                                                                                      processed_rds,
                                                                                                                                                      singler_model_file
                                                                                                                                                      )}

      // creates input for cellassign [meta, cellassign ref file, cell assign ref name]
      cellassign_input_ch = grouped_celltype_ch
        .map{meta, processed_hdf5, unfiltered_rds, filtered_rds, processed_rds, singler_model_file, cellassign_ref_file, cellassign_ref_name -> tuple(meta,
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

      // get CellAssign cell type predictions and add them to SCE
      classify_cellassign(all_celltype_assignments_ch)

      emit: classify_cellassign.out

}
