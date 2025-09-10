// run inferCNV on an SCE object that has consensus cell types

process call_infercnv {
  container params.SCPCATOOLS_INFERCNV_CONTAINER
  publishDir (
    path: "${meta.infercnv_dir}",
    mode: 'copy',
    pattern: "${meta.unique_id}_infercnv-*"
  )
  label 'mem_32'
  label 'cpus_8'
  tag "${meta.unique_id}"
  input:
    tuple val(meta), path(processed_rds), path(infercnv_gene_order)
  output:
    tuple val(meta), path(processed_rds), path(results_file), path(heatmap_file)
  script:
    results_file="${meta.unique_id}_infercnv-results.rds"
    heatmap_file="${meta.unique_id}_infercnv-heatmap.png"
    """
    touch ${heatmap_file}
    touch ${results_file}
    #run_infercnv.R \
    #  --input_sce_file ${processed_rds} \
    #  --output_rds ${results_file} \
    #  --output_heatmap ${heatmap_file} \
    #  --temp_dir \$PWD \
    #  --gene_order_file ${infercnv_gene_order} \
    #  --threads ${task.cpus} \
    #  ${params.seed ? "--random_seed ${params.seed}" : ""}
    """
  stub:
    results_file="${meta.unique_id}_infercnv-results.rds"
    heatmap_file="${meta.unique_id}_infercnv-heatmap.png"
    """
    touch "${results_file}"
    touch "${heatmap_file}"
    """
}



process add_infercnv_to_sce {
  container params.SCPCATOOLS_SLIM_CONTAINER
  label 'mem_8'
  tag "${meta.library_id}"
  input:
    tuple val(meta), path(processed_rds), path(infercnv_results_file)
  output:
    tuple val(meta), path(infercnv_sce)
  script:
    // call this sce to avoid confusing with the infercnv results rds file
    infercnv_sce = "${processed_rds.baseName}_infercnv.rds"
    """
    cp ${processed_rds} ${infercnv_sce}

    #add_infercnv_to_sce.R \
    #  --input_sce_file ${processed_rds} \
    #  --infercnv_results_file "${infercnv_results_file}" \
    #  --output_sce_file ${infercnv_sce}
    """
  stub:
    infercnv_sce = "${processed_rds.baseName}_infercnv.rds"
    """
    touch "${infercnv_sce}"
    """
}


workflow run_infercnv {
  take: sce_files_channel // channel of meta, unfiltered_sce, filtered_sce, processed_sce
  main:
    def empty_file = "${projectDir}/assets/NO_FILE"
    // read in sample metadata and make a list of cell line samples; we won't run inferCNV on these
    cell_line_samples = Channel.fromPath(params.sample_metafile)
      .splitCsv(header: true, sep: '\t')
      .map{
        [
          sample_id: it.scpca_sample_id,
          is_cell_line: Utils.parseNA(it.is_cell_line).toBoolean()
        ]
      }
      .filter{it.is_cell_line}
      .map{it.sample_id}
      .toList()

    sce_files_channel_branched = sce_files_channel
      .branch{
          not_eligible: (
            it[0]["sample_id"].split(",").collect{it in cell_line_samples.getVal()}.every()
            || it[0]["infercnv_reference_cell_count"] < params.infercnv_min_reference_cells
          )
          eligible: true
      }

    // create input for infercnv: [augmented meta, processed_sce, gene order file]
    infercnv_prepared_ch = sce_files_channel_branched.eligible
      .map{ meta_in, unfiltered_sce, filtered_sce, processed_sce ->
        def meta = meta_in.clone(); // local copy for safe modification
        meta.infercnv_dir = "${params.checkpoints_dir}/infercnv/${meta.unique_id}";
        meta.infercnv_heatmap_file = "${meta.infercnv_dir}/${meta.unique_id}_infercnv-heatmap.png";
        meta.infercnv_results_file = "${meta.infercnv_dir}/${meta.unique_id}_infercnv-results.rds";
        // return simplified input with gene order file
        [meta, processed_sce, file("${meta.infercnv_gene_order}", checkIfExists: true)]
      }

    // branch for run conditions
    infercnv_input_ch = infercnv_prepared_ch
      .branch{
        skip_infercnv: (
          !params.repeat_infercnv
          && file(it[0].infercnv_heatmap_file).exists()
          && file(it[0].infercnv_results_file).exists()
        )
        run_infercnv: true
      }

    // run inferCNV
    // outputs: [meta, processed sce, results file, heatmap]
    call_infercnv(infercnv_input_ch.run_infercnv)

    // prepare to add results for all eligible libraries: [meta, processed sce, results file]
    add_infercnv_results_ch = infercnv_input_ch.skip_infercnv
      .map{ meta, processed_sce, gene_order_file -> tuple(
        meta,
        processed_sce,
        file("${meta.infercnv_results_file}", checkIfExists: true),
        file("${meta.infercnv_heatmap_file}", checkIfExists: true)
      )}
      .mix(call_infercnv.out)
      .map{[it[0], it[1], it[2]]} // remove heatmap

    add_infercnv_to_sce(add_infercnv_results_ch)

    // add back in the unchanged sce files to the results
    export_channel = add_infercnv_to_sce.out
      .map{meta, processed_sce -> tuple(
        meta.unique_id,
        meta,
        processed_sce
        )}
      // add in unfiltered and filtered sce files, for eligible samples only
      .join(
        sce_files_channel_branched.eligible.map{[it[0]["unique_id"], it[1], it[2]]},
        by: 0, failOnMismatch: true, failOnDuplicate: true
      )
      // rearrange back to [meta, unfiltered, filtered, processed]
      .map{_unique_id, meta, processed_sce, unfiltered_sce, filtered_sce ->
        [meta, unfiltered_sce, filtered_sce, processed_sce]
      }
      // mix in libraries which we did not run inferCNV on
      .mix(sce_files_channel_branched.not_eligible)

    emit: export_channel

  }
