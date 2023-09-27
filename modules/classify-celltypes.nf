
process classify_singleR {
    container params.SCPCATOOLS_CONTAINER
    publishDir (
        path: "${params.checkpoints_dir}/celltype/${meta.library_id}",
        mode:  'copy',
        pattern: "*{results.rds,.tsv,.json}" // Everything except processed rds
    )
    label 'mem_8'
    label 'cpus_4'
    input:
        tuple val(meta), path(processed_rds), path(singler_model_file)
    output:
        tuple val(meta), path(processed_rds), path(singler_annotations_tsv), path(singler_full_results)
    script:

      singler_annotations_tsv = "${meta.library_id}_singler_annotations.tsv"
      singler_full_results = "${meta.library_id}_singler_full_results.rds"
      """
      classify_SingleR.R \
        --sce_file ${processed_rds} \
        --singler_model_file ${singler_model_file} \ # todo: filename should contain celldex version
        --output_singler_annotations_file ${singler_annotations_tsv} \
        --output_singler_results_file ${singler_full_results} \
        --seed ${params.seed} \
        --threads ${task.cpus}
      """
    stub:
      singler_annotations_tsv = "${meta.library_id}_singler_annotations.tsv"
      singler_full_results = "${meta.library_id}_singler_full_results.rds"
      """
      touch "${singler_annotations_tsv}"
      touch "${singler_full_results}"
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
          // add ref name for cellassign since we cannot store it in the cellassign output
          // singler ref name does not need to be added because it is stored in the singler model
          cellassign_ref_name = it.cellassign_ref_name
        ]}

      // create input for singleR: [meta, processed, SingleR reference model]
      singler_input_ch = processed_sce_channel
        .map{[it[0]["project_id"]] + it}
        .combine(celltype_ch, by: 0)
        .map{it.drop(1)} // remove extra project ID

      // perform singleR celltyping and export TSV
      classify_singleR(singler_input_ch)

    // temporary during development
    emit: classify_singleR.out

}
