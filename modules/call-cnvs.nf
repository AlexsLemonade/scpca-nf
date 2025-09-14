// run inferCNV on an SCE object that has consensus cell types

process run_infercnv {
  container params.SCPCATOOLS_INFERCNV_CONTAINER
  publishDir (
    path: "${meta.infercnv_dir}",
    mode: 'copy',
    pattern: "${meta.unique_id}_infercnv-*"
  )
  label 'mem_96'
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
    # If we don't have sufficient cells to run inferCNV, just touch the files
    # must use < to accommodate a possible null value; null < X is always true
    if [[ ${meta.infercnv_reference_cell_count} -lt ${params.infercnv_min_reference_cells} ]]; then
      touch "${results_file}"
      touch "${heatmap_file}"
    else
      # note that inferCNV fails, the script will output empty results/heatmap files
      run_infercnv.R \
        --input_sce_file ${processed_rds} \
        --output_rds ${results_file} \
        --output_heatmap ${heatmap_file} \
        --temp_dir \$PWD \
        --gene_order_file ${infercnv_gene_order} \
        --threads ${task.cpus} \
        ${params.seed ? "--random_seed ${params.seed}" : ""}
    fi
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
  tag "${meta.unique_id}"
  input:
    tuple val(meta), path(processed_rds), path(infercnv_results_file)
  output:
    tuple val(meta.unique_id), val(meta), path(infercnv_sce)
  script:
    // call this sce to avoid confusing it with the infercnv_results_file rds
    infercnv_sce = "${processed_rds.baseName}_infercnv.rds"
    """
    add_infercnv_to_sce.R \
      --input_sce_file ${processed_rds} \
      --infercnv_results_file "${infercnv_results_file}" \
      --infercnv_threshold ${params.infercnv_min_reference_cells} \
      --output_sce_file ${infercnv_sce}
    """
  stub:
    infercnv_sce = "${processed_rds.baseName}_infercnv.rds"
    """
    touch "${infercnv_sce}"
    """
}


workflow call_cnvs {
  take: sce_files_channel // channel of meta, unfiltered_sce, filtered_sce, processed_sce
  main:
    // read in sample metadata and make a list of cell line samples; we won't run inferCNV on these
    cell_line_samples = Channel.fromPath(params.sample_metafile)
      .splitCsv(header: true, sep: '\t')
      .filter{it.is_cell_line.toBoolean()}
      .map{it.scpca_sample_id}
      .toList()

    sce_files_channel_branched = sce_files_channel
      .branch{
          cell_line: it[0]["sample_id"].split(",").every{it in cell_line_samples.getVal()}
          tissue: true
      }

    // create input for infercnv: [augmented meta, processed_sce, gene order file]
    infercnv_prepared_ch = sce_files_channel_branched.tissue
      .map{ meta_in, unfiltered_sce, filtered_sce, processed_sce ->
        def meta = meta_in.clone(); // local copy for safe modification
        meta.infercnv_dir = "${params.checkpoints_dir}/infercnv/${meta.unique_id}";
        meta.infercnv_heatmap_file = "${meta.infercnv_dir}/${meta.unique_id}_infercnv-heatmap.png";
        meta.infercnv_results_file = "${meta.infercnv_dir}/${meta.unique_id}_infercnv-results.rds";
        // return simplified input with gene order file
        [meta, processed_sce, file("${meta.infercnv_gene_order}", checkIfExists: true)]
      }

    // branch to skip libraries with existing results if repeat is off
    infercnv_input_ch = infercnv_prepared_ch
      .branch{
        skip_infercnv: (
          !params.repeat_cnv_inference
          && file(it[0].infercnv_heatmap_file).exists()
          && file(it[0].infercnv_results_file).exists()
        )
        call_infercnv: true
      }

    infercnv_input_ch.call_infercnv.view()

    // run inferCNV
    // outputs: [meta, processed sce, results file, heatmap]
    run_infercnv(infercnv_input_ch.call_infercnv)

    // prepare to add results for all eligible libraries: [meta, processed sce, results file]
    add_infercnv_results_ch = infercnv_input_ch.skip_infercnv
      .map{ meta, processed_sce, gene_order_file -> tuple(
        meta,
        processed_sce,
        file("${meta.infercnv_results_file}", checkIfExists: true),
        file("${meta.infercnv_heatmap_file}", checkIfExists: true)
      )}
      .mix(run_infercnv.out)
      // drop the heatmap file
      .map{it.dropRight(1)}

    // add inferCNV results to the SCE object
    // returns [ unique id, meta, processed sce ]
    add_infercnv_to_sce(add_infercnv_results_ch)

    export_channel = add_infercnv_to_sce.out
      // add in unfiltered and filtered sce files, for eligible samples only
      .join(
        sce_files_channel_branched.tissue.map{[it[0]["unique_id"], it[1], it[2]]},
        by: 0, failOnMismatch: true, failOnDuplicate: true
      )
      // rearrange back to [meta, unfiltered, filtered, processed]
      .map{_unique_id, meta, processed_sce, unfiltered_sce, filtered_sce ->
        [meta, unfiltered_sce, filtered_sce, processed_sce]
      }
      // mix in libraries which we did not run inferCNV on
      .mix(sce_files_channel_branched.cell_line)

    emit: export_channel

  }
