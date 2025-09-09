// run inferCNV on an SCE object that has consensus cell types
process call_infercnv {
  container params.SCPCATOOLS_INFERCNV_CONTAINER
  label 'mem_16'
  tag "${meta.library_id}"
  input:
    tuple val(meta), path(unfiltered_rds), path(filtered_rds), path(processed_rds), path(infercnv_gene_order)
  output:
    tuple val(meta), path(unfiltered_rds), path(filtered_rds), path(infercnv_rds)
  script:
    infercnv_dir = file(meta.infercnv_checkpoints_dir).name
    infercnv_rds = "${processed_rds.baseName}_infercnv.rds"
    """
    # create output directory to store the infercnv.png file
    mkdir -p ${infercnv_dir}

    run_infercnv.R \
      --input_sce_file ${processed_rds} \
      --output_sce_file ${infercnv_rds} \
      --output_dir ${infercnv_dir} \
      --gene_order_file ${infercnv_gene_order} \
      --num_threads ${task.cpus} \
      ${params.seed ? "--random_seed ${params.seed}" : ""}
    """
  stub:
    infercnv_dir = file(meta.infercnv_checkpoints_dir).name
    infercnv_rds = "${processed_rds.baseName}_infercnv.rds"
    """
    mkdir "${infercnv_dir}"
    touch "${infercnv_dir}/infercnv.png"
    touch "${infercnv_rds}"
    """
}



workflow run_infercnv {
  take: sce_files_channel // channel of meta, unfiltered_sce, filtered_sce, processed_sce
  main:
    // get cell line samples, which we will not run infercnv on since they do not get cell typed
    cell_line_samples = Channel.fromPath(params.sample_metafile)
      .splitCsv(header: true, sep: '\t')
      .map{
        [
          sample_id: it.scpca_sample_id,
          is_cell_line: Utils.parseNA(it.is_cell_line).toBoolean() // FALSE -> false, NA -> false, TRUE -> true
        ]
      }
      .filter{it.is_cell_line}
      .map{it.sample_id}
      .toList()

    sce_files_channel_branched = sce_files_channel
      // add in infercnv output directory and png heatmap file for later use
      // do this _before_ branching so the heatmap file can always be checked when prepping the report
      .map{ meta_in, unfiltered_sce, filtered_sce, processed_sce ->
          def meta = meta_in.clone(); // local copy for safe modification
          meta.infercnv_checkpoints_dir = "${params.checkpoints_dir}/infercnv/${meta.library_id}";
          meta.infercnv_heatmap_file = "${meta.infercnv_checkpoints_dir}/infercnv.png"; // using inferCNV's naming
          [meta, unfiltered_sce, filtered_sce, processed_sce]
      }
      // branch to skip if:
      //  - repeat_infercnv is false and infercnv results already exist
      //  - cell lines
      //  - libraries where number of cells is below the threshold
      .branch{
        skip_infercnv: (
          (!params.repeat_infercnv && file(it[0].infercnv_heatmap_file).exists())
          || it[0]["sample_id"].split(",").collect{it in cell_line_samples.getVal()}.every()
          || it[0]["infercnv_reference_cell_count"] < params.infercnv_min_reference_cells
        )
        run_infercnv: true
      }

    infercnv_ch = sce_files_channel_branched.run_infercnv
      // add in the gene order file
      .map{it.toList() + [file(it[0].infercnv_gene_order)]}

    // run inferCNV
    call_infercnv(infercnv_ch)

    // mix with skipped libraries
    combined_ch = call_infercnv.out
      .mix(sce_files_channel_branched.skip_infercnv)

    emit: combined_ch
}
