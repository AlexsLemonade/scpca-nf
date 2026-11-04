// run inferCNV on an SCE object that has consensus cell types

include { getVersions; makeJson; getMetaVal; checkpointIsReusable; pullthroughContainer } from '../lib/utils.nf'

process run_infercnv {
  container "${pullthroughContainer(params.scpcatools_infercnv_container, params.pullthrough_registry)}"
  publishDir (
    path: "${meta.infercnv_checkpoints_dir}",
    mode: 'copy',
    pattern: "{${meta.unique_id}_infercnv-*,scpca-meta.json}"
  )
  label 'mem_128'
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

    meta += getVersions()
    meta_json = makeJson(meta)
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
    meta += getVersions()
    meta_json = makeJson(meta)
    """
    touch "${results_file}"
    touch "${table_file}"
    touch "${heatmap_file}"
    echo '${meta_json}' > scpca-meta.json
    """
}


process add_infercnv_to_sce {
  container "${pullthroughContainer(params.scpcatools_slim_container, params.pullthrough_registry)}"
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
    cell_line_samples = channel.fromPath(params.sample_metafile)
      .splitCsv(header: true, sep: '\t')
      .filter{ it.is_cell_line.toBoolean() }
      .map{ it -> it.scpca_sample_id }
      .toList()

    // branch to run inferCNV on tissue samples with non-null reference cell counts only
    sce_files_channel_branched = sce_files_channel
      .branch{ meta, _unfiltered, _filtered, _processed ->
        no_infercnv: (
          meta.sample_id.split(",").every{ it in cell_line_samples.getVal() }
          || meta.infercnv_reference_cell_count == null
        )
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
        def gene_order_file = meta.infercnv_gene_order ? file(meta.infercnv_gene_order, checkIfExists: true) : []

        if (!gene_order_file) {
          log.warn("Gene order file for '${meta.unique_id}' was not found. inferCNV will be skipped for this library.")
        }

        // return simplified input with gene order file
        [meta, processed_sce, gene_order_file]
      }


    // branch libraries into one of three paths:
    // - cannot_run: inferCNV prerequisites are not met (not enough normal reference cells or no gene order file)
    //   these libraries get empty placeholder results, and any existing checkpoint files are ignored
    // - reuse_checkpoint: repeat is off and there are unchanged existing results
    //   (if this is a release run, existing results must also have been created by a release)
    // - run_infercnv: everything else
    infercnv_input_ch = infercnv_prepared_ch
      .branch{ meta, _processed_sce, gene_order_file ->
        def stored_cell_hash = getMetaVal(file("${meta.infercnv_checkpoints_dir}/scpca-meta.json"), "infercnv_reference_cell_hash")

        cannot_run: (
          meta.infercnv_reference_cell_count < params.infercnv_min_reference_cells
          || !(gene_order_file in Path && gene_order_file.exists()) // if the gene order file doesn't exist, we can't run infercnv
        )
        reuse_checkpoint: (
          !params.repeat_cnv_inference
          && file(meta.infercnv_heatmap_file).exists()
          && file(meta.infercnv_results_file).exists()
          && file(meta.infercnv_table_file).exists()
          && meta.infercnv_reference_cell_hash == stored_cell_hash
          && checkpointIsReusable(file("${meta.infercnv_checkpoints_dir}/scpca-meta.json"))
        )
        run_infercnv: true
      }

    // run inferCNV
    // outputs: [meta, processed sce, results file, table file, heatmap file]
    run_infercnv(infercnv_input_ch.run_infercnv)

    // existing checkpoint results: [meta, processed sce, results file, table file, heatmap file]
    reuse_infercnv_ch = infercnv_input_ch.reuse_checkpoint
      .map{ meta, processed_sce, _gene_order_file ->
        [meta, processed_sce, file(meta.infercnv_results_file), file(meta.infercnv_table_file), file(meta.infercnv_heatmap_file)]
      }

    // empty placeholder results for libraries where inferCNV cannot be run
    // placeholders are written to the work directory so existing checkpoint files are never used or overwritten
    // outputs: [meta, processed sce, results file, table file, heatmap file]
    empty_infercnv_ch = infercnv_input_ch.cannot_run
      .map{ meta, processed_sce, _gene_order_file ->
        def placeholder_dir = file("${workflow.workDir}/infercnv-placeholders/${meta.unique_id}")
        placeholder_dir.mkdirs()
        def infercnv_results = placeholder_dir / "${meta.unique_id}_infercnv-results.rds"
        def infercnv_table   = placeholder_dir / "${meta.unique_id}_infercnv-table.txt"
        def infercnv_heatmap = placeholder_dir / "${meta.unique_id}_infercnv-heatmap.png"
        // always write empty files, so add_infercnv_to_sce.R sets the status based on the reference cell count
        // https://docs.seqera.io/nextflow/working-with-files#reading-and-writing-an-entire-file
        infercnv_results.text = ''
        infercnv_table.text = ''
        infercnv_heatmap.text = ''

        [meta, processed_sce, infercnv_results, infercnv_table, infercnv_heatmap]
      }

    // combine all libraries to add inferCNV results to the SCE
    add_infercnv_results_ch = reuse_infercnv_ch
      .mix(empty_infercnv_ch, run_infercnv.out.infercnv)


    // add inferCNV results to the SCE object
    // returns [ meta, processed sce, infercnv_heatmap_file ]
    // note we keep the heatmap file to eventually stage for QC report
    add_infercnv_to_sce(add_infercnv_results_ch)

    export_channel = add_infercnv_to_sce.out
      .map{ meta, processed_sce, infercnv_heatmap_file ->
        [meta.unique_id, meta, processed_sce, infercnv_heatmap_file]
      }
      // add in unfiltered and filtered sce files, for tissue samples only
      .join(
        sce_files_channel_branched.tissue.map{ meta, unfiltered, filtered, _processed ->
          [meta.unique_id, unfiltered, filtered]
        },
        by: 0, failOnMismatch: true, failOnDuplicate: true
      )
      // rearrange back to [meta, unfiltered, filtered, processed, infercnv_heatmap_file]
      .map{ _unique_id, meta, processed_sce, heatmap_file, unfiltered_sce, filtered_sce ->
        [meta, unfiltered_sce, filtered_sce, processed_sce, heatmap_file]
      }
      // mix in libraries which we did not run inferCNV on
      .mix(
        // add in an empty file for heatmap placeholder first
        sce_files_channel_branched.no_infercnv
          .map{ meta, unfiltered_sce, filtered_sce, processed_sce ->
            [meta, unfiltered_sce, filtered_sce, processed_sce, []]
          }
      )

    emit: export_channel

  }
