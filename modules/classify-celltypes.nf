
process classify_singleR {
    container params.SCPCATOOLS_CONTAINER
    publishDir (
        path: "${params.checkpoints_dir}/celltype/${meta.library_id}",
        mode:  'copy',
        pattern: "*{_singler.rds,.tsv,.json}" // Everything except processed rds
    )
    label 'mem_8'
    label 'cpus_4'
    input:
        tuple val(meta), path(processed_rds), path(singler_model_file)
    output:
        tuple val(meta), path(processed_rds), path(singler_annotations_tsv), path(singler_full_results)
    script:
      singler_annotations_tsv = "${meta.library_id}_singler_annotations.tsv"
      singler_full_results = "${meta.library_id}_singler.rds"
      """
      classify_SingleR.R \
        --sce_file ${processed_rds} \
        --singler_model_file ${singler_model_file} \
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


process classify_cellassign {
  container params.SCPCATOOLS_CONTAINER
    publishDir (
        path: "${params.checkpoints_dir}/celltype/${meta.library_id}",
        mode:  'copy',
        pattern: "*{.tsv,.json}" // Only the prediction matrix (tsv) and meta
    )
  label 'mem_32'
  label 'cpus_12'
  input:
    tuple val(meta), path(processed_rds), path(cellassign_ref)
  output:
    tuple val(meta), path(processed_rds), path(cellassign_predictions_tsv)
  script:
    processed_hdf5 = "${meta.library_id}_processed.hdf5"
    cellassign_predictions = "${meta.library_id}_predictions.tsv"

    """
    # Convert SCE to AnnData
    sce_to_anndata.R \
        --input_sce_file ${processed_rds} \
        --output_rna_h5 ${processed_hdf5} 
    
    # Run CellAssign
    predict_cellassign.py \
      --input_hdf5_file ${processed_hdf5} \
      --output_predictions ${cellassign_predictions_tsv} \
      --reference ${cellassign_ref} \
      --seed ${params.seed} \
      --threads ${task.cpus}
    """
  stub:
    processed_hdf5 = "${meta.library_id}_processed.hdf5"
    cellassign_predictions = "${meta.library_id}_predictions.tsv"
    """
    touch "${cellassign_predictions}"
    touch "${processed_hdf5}"
    """
}

// TODO - overhaul this process
process add_celltypes {
  container params.SCPCATOOLS_CONTAINER
  publishDir "${params.results_dir}/${meta.project_id}/${meta.sample_id}", mode: 'copy'
  label 'mem_4'
  label 'cpus_2'
  input:
    tuple val(meta), path(input_rds), path(cellassign_predictions), val(ref_name)
  output:
    tuple val(meta), path(annotated_rds)
  script:
    annotated_rds = "${meta.library_id}_annotated.rds"
    """
    classify_cellassign.R \
      --input_sce_file ${input_rds} \
      --output_sce_file ${annotated_rds} \
      --cellassign_predictions ${cellassign_predictions} \
      --reference_name ${ref_name}
    """
}

workflow annotate_celltypes {
    take: sce_files_channel // channel of meta, unfiltered_sce, filtered_sce, processed_sce
    main:
      // get just the meta and processed sce
      processed_sce_channel = sce_files_channel.map{[it[0], it[3]]}

      // channel with celltype model and project ids
      celltype_ch = Channel.fromPath(params.celltype_project_metafile)
        .splitCsv(header: true, sep: '\t')
        .map{[
         // project id
         it.scpca_project_id,
         // singler model file
         Utils.parseNA(it.singler_ref_file) ? file("${params.singler_models_dir}/${it.singler_ref_file}") : null,
         // cellassign reference file
         Utils.parseNA(it.cellassign_ref_file) ? file("${params.cellassign_ref_dir}/${it.cellassign_ref_file}") : null,
         // add ref name for cellassign since we cannot store it in the cellassign output
         // singler ref name does not need to be added because it is stored in the singler model
         Utils.parseNA(it.cellassign_ref_name)
        ]}

      // create cell typing channel: [meta, processed_rds, singler_model, cellassign_ref, cellassign_ref_name]
      celltype_input_ch = processed_sce_channel
        .map{[it[0]["project_id"]] + it}
        .combine(celltype_ch, by: 0)
        .map{it.drop(1)} // remove extra project ID
        // we only run celltyping for rows with a singler model file
        // branch here so we have meta and processed sce in the .skip
        .branch{
          skip: it[2] == null
          run: true
        }
        
      
      // create input for singleR: [meta, processed, SingleR reference model]
      singler_input_ch = celltype_input_ch.run
        .map{meta, processed_rds, singler_model, cellassign_ref, cellassign_ref_name ->
             [meta, processed_rds, singler_model]}
             
      // perform singleR celltyping and export TSV
      classify_singleR(singler_input_ch) 

      // proceed to cellassign, only for datasets with references
      cellassign_ch = classify_singleR.out //  meta, processed_rds, singler_annotations_tsv singler_full_results
        .join(celltype_input_ch.run.map{[it[0], it[3]]}, // meta, cellassign_ref
              by: 0, failOnMismatch: true, failOnDuplicate: true)
       .map{meta, processed_rds, singler_annotations_tsv, singler_full_results, cellassign_ref ->
          [meta, processed_rds, cellassign_ref]}
      
     // classify_cellassign(cellassign_ch) // output: meta, processed_rds, cellassign_predictions, cellassign_ref_name

      // TODO mix celltyping results back up with `celltype_input_ch.skip`

      // add back in the unchanged sce files
      // TODO update below with output channel results:
      export_channel = processed_sce_channel
        .map{[it[0]["library_id"]] + it}
        // add in unfiltered and filtered sce files
        .join(sce_files_channel.map{[it[0]["library_id"], it[1], it[2]]},
              by: 0, failOnMismatch: true, failOnDuplicate: true)
        // rearrange to be [meta, unfiltered, filtered, processed]
        .map{library_id, meta, processed_sce, unfiltered_sce, filtered_sce ->
            [meta, unfiltered_sce, filtered_sce, processed_sce]}

    emit: export_channel

}
