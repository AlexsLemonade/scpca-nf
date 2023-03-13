
process classify_singleR {
    container params.SCPCATOOLS_CONTAINER
    label 'mem_8'
    label 'cpus_4'
    publishDir "${params.results_dir}/${meta.project_id}/${meta.sample_id}", mode: 'copy'
    input:
        tuple val(meta), path(processed_rds), path(singler_models)
    output:
        tuple val(meta), path(annotated_rds)
    script:
      annotated_rds = "${meta.library_id}_annotated.rds"
      singler_list = singler_models in List ? singler_models.join(',') : singler_models
      """
      classify_SingleR.R \
        --input_sce_file ${processed_rds} \
        --output_sce_file ${annotated_rds} \
        --singler_models ${singler_list} \
        --seed ${params.seed} \
        --threads ${task.cpus}

      """


}

workflow annotate_celltypes {
    take: processed_sce_channel
    main:
      // channel with celltype model and project ids
      celltype_ch = Channel.fromPath(params.project_celltype_metafile)
        .splitCsv(header: true, sep: '\t')
        .map{[
          project_id = it.scpca_project_id,
          singler_model = "${params.celltype_model_dir}/${it.celltype_ref_name}_model.rds"
        ]}

      // create channel with grouped meta, processed sce object, and all references to use
      grouped_celltype_ch = processed_sce_channel
        .map{[it[0]["project_id"]] + it}
        .combine(celltype_ch, by: 0)
        .map{it.drop(1)} // remove extra project ID
        .groupTuple(by: 0) // group by meta
        .map{[
          it[0], // meta
          it[1][0], // processed rds
          it[2] // tuple of reference models
        ]}

      classify_singleR(grouped_celltype_ch)

    emit: classify_singleR.out

}
