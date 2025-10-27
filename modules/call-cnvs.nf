// run inferCNV on an SCE object that has consensus cell types

process run_infercnv {
  container params.SCPCATOOLS_INFERCNV_CONTAINER
  publishDir (
    path: "${meta.infercnv_checkpoints_dir}",
    mode: 'copy',
    pattern: "{${meta.unique_id}_infercnv-*,scpca-meta.json}"
  )
  label 'mem_96'
  label 'cpus_8'
  tag "${meta.unique_id}"
  input:
    tuple val(meta), path(processed_rds), path(infercnv_gene_order)
  output:
    tuple val(meta), path(processed_rds), path(results_file), path(table_file), path(heatmap_file), emit: infercnv
    path "scpca-meta.json", emit: metafile
  script:
    results_file="${meta.unique_id}_infercnv-results.rds"
    table_file="${meta.unique_id}_infercnv-table.txt"
    heatmap_file="${meta.unique_id}_infercnv-heatmap.png"

    meta_json = Utils.makeJson(meta)
    """
    # note that if inferCNV fails, the script will output empty results/heatmap files
    mkdir infercnv_tmp
    run_infercnv.R \
      --input_sce_file ${processed_rds} \
      --output_rds ${results_file} \
      --output_table ${table_file} \
      --output_heatmap ${heatmap_file} \
      --temp_dir "infercnv_tmp" \
      --gene_order_file ${infercnv_gene_order} \
      --threads ${task.cpus} \
      ${params.seed ? "--random_seed ${params.seed}" : ""}

    # write out meta file
    echo '${meta_json}' > scpca-meta.json
    """
  stub:
    results_file="${meta.unique_id}_infercnv-results.rds"
    table_file="${meta.unique_id}_infercnv-table.txt"
    heatmap_file="${meta.unique_id}_infercnv-heatmap.png"
    meta_json = Utils.makeJson(meta)
    """
    touch "${results_file}"
    touch "${table_file}"
    touch "${heatmap_file}"
    echo '${meta_json}' > scpca-meta.json
    """
}


process add_infercnv_to_sce {
  container params.SCPCATOOLS_SLIM_CONTAINER
  label 'mem_8'
  tag "${meta.unique_id}"
  input:
    tuple val(meta), path(processed_rds), path(infercnv_results_file), path(infercnv_table_file), path(infercnv_heatmap_file)
  output:
    tuple val(meta), path(infercnv_sce), path(infercnv_heatmap_file)
  script:
    // call this sce to avoid confusing it with the infercnv rds results file
    infercnv_sce = "${processed_rds.baseName}_infercnv.rds"
    """
    add_infercnv_to_sce.R \
      --input_sce_file ${processed_rds} \
      --infercnv_results_file "${infercnv_results_file}" \
      --infercnv_table_file "${infercnv_table_file}" \
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
      .filter{ it.is_cell_line.toBoolean() }
      .map{ it -> it.scpca_sample_id }
      .toList()

    sce_files_channel_branched = sce_files_channel
      .branch{ it ->
        cell_line: it[0]["sample_id"].split(",").every{ samples -> samples[0] in cell_line_samples.getVal() }
        tissue: true
      }

    // create input for infercnv: [augmented meta, processed_sce, gene order file]
    infercnv_prepared_ch = sce_files_channel_branched.tissue
      .map{ meta_in, _unfiltered_sce, _filtered_sce, processed_sce ->
        def meta = meta_in.clone() // local copy for safe modification
        // define infercnv checkpoint files
        meta.infercnv_checkpoints_dir = "${params.checkpoints_dir}/infercnv/${meta.unique_id}"
        meta.infercnv_heatmap_file = "${meta.infercnv_checkpoints_dir}/${meta.unique_id}_infercnv-heatmap.png"
        meta.infercnv_results_file = "${meta.infercnv_checkpoints_dir}/${meta.unique_id}_infercnv-results.rds"
        meta.infercnv_table_file = "${meta.infercnv_checkpoints_dir}/${meta.unique_id}_infercnv-table.txt"
        // if meta.infercnv_reference_cell_count doesn't exist, set it to null
        // this happens when perform_celltyping is on but a library has no cell type references
        meta.infercnv_reference_cell_count = meta.infercnv_reference_cell_count ?: null

        // return simplified input with gene order file
        [meta, processed_sce, file(meta.infercnv_gene_order, checkIfExists: true)]
      }


    // branch to skip libraries when either:
    // - repeat is off and there are unchanged existing results
    // - there are not enough normal reference cells
    infercnv_input_ch = infercnv_prepared_ch
      .branch{ it ->
        def stored_cell_hash = Utils.getMetaVal(file("${it[0].infercnv_checkpoints_dir}/scpca-meta.json"), "infercnv_reference_cell_hash")

        skip_infercnv: (
        (
          !params.repeat_cnv_inference
          && file(it[0].infercnv_heatmap_file).exists()
          && file(it[0].infercnv_results_file).exists()
          && file(it[0].infercnv_table_file).exists()
          && it[0].infercnv_reference_cell_hash == stored_cell_hash
        ) || it[0].infercnv_reference_cell_count < params.infercnv_min_reference_cells
        )
        run_infercnv: true
      }

    // run inferCNV
    // outputs: [meta, processed sce, results file, table file, heatmap file]
    run_infercnv(infercnv_input_ch.run_infercnv)

    // prepare to add results for all eligible libraries: [meta, processed sce, results file]
    add_infercnv_results_ch = infercnv_input_ch.skip_infercnv
      .map{ meta, processed_sce, _gene_order_file ->
        def infercnv_results = file(meta.infercnv_results_file)
        def infercnv_table   = file(meta.infercnv_table_file)
        def infercnv_heatmap = file(meta.infercnv_heatmap_file)
        // ensure the infercnv_checkpoints_dir exists before proceeding; all these files are in the same dir
        infercnv_results.parent.mkdirs()
        // create empty files if they don't exist
        // https://www.nextflow.io/docs/latest/working-with-files.html#reading-and-writing-an-entire-file
        if (!infercnv_results.exists()) infercnv_results.text = ''
        if (!infercnv_table.exists())   infercnv_table.text = ''
        if (!infercnv_heatmap.exists()) infercnv_heatmap.text = ''

        // return simplified input with gene order file
        [meta, processed_sce, infercnv_results, infercnv_table, infercnv_heatmap]
      }
      .mix(run_infercnv.out.infercnv)


    // add inferCNV results to the SCE object
    // returns [ meta, processed sce, infercnv_heatmap_file ]
    // note we keep the heatmap file to eventually stage for QC report
    add_infercnv_to_sce(add_infercnv_results_ch)


    export_channel = add_infercnv_to_sce.out
      .map{ meta, processed_sce, infercnv_heatmap_file ->
        [meta["unique_id"], meta, processed_sce, infercnv_heatmap_file]
      }
      // add in unfiltered and filtered sce files, for tissue samples only
      .join(
        sce_files_channel_branched.tissue.map{ meta, unfiltered, filtered, _processed ->
          [meta["unique_id"], unfiltered, filtered]
        },
        by: 0, failOnMismatch: true, failOnDuplicate: true
      )
      // rearrange back to [meta, unfiltered, filtered, processed, infercnv_heatmap_file]
      .map{ it -> it.drop(1) }
      // mix in cell line libraries which we did not run inferCNV on
      .mix(
        // add in an empty file for heatmap placeholder first
        sce_files_channel_branched.cell_line
          .map{ meta, unfiltered_sce, filtered_sce, processed_sce ->
            [meta, unfiltered_sce, filtered_sce, processed_sce, []]
          }
      )

    emit: export_channel

  }
