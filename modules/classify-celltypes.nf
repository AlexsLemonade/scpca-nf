
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

workflow annotate_celltypes {
    take: processed_sce_channel
    main:
      // channel with celltype model and project ids
      celltype_ch = Channel.fromPath(params.celltype_refs_metafile)
        .splitCsv(header: true, sep: '\t')
        .map{[
          project_id = it.scpca_project_id,
          singler_model_file = "${params.celltype_model_dir}/${it.celltype_ref_name}_model.rds"
        ]}

      // create channel grouped_celltype_ch as: [meta, processed sce object, SingleR reference model]
      // input processed_sce_channel is [meta, unfiltered, filtered, processed]
      grouped_celltype_ch = processed_sce_channel
        .map{[it[0]["project_id"]] + it}
        .combine(celltype_ch, by: 0)
        .map{it.drop(1)} // remove extra project ID
        // creates [meta, processed, SingleR reference model]
        .map{ it[0..2] }

      classify_singleR(grouped_celltype_ch)

    emit: classify_singleR.out

}
