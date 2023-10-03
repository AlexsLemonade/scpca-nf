
process classify_singler {
    container params.SCPCATOOLS_CONTAINER
    publishDir (
      path: "${meta.celltype_publish_dir}",
      mode: 'copy',
      pattern: "${singler_dir}"
    )
    label 'mem_8'
    label 'cpus_4'
    input:
      tuple val(meta), path(processed_rds), path(singler_model_file)
    output:
      tuple val(meta.library_id), path(singler_dir)
    script:
      singler_dir = file(meta.singler_dir).name
      """
      # create output directory
      mkdir "${singler_dir}"

      classify_SingleR.R \
        --sce_file "${processed_rds}" \
        --singler_model_file "${singler_model_file}" \
        --output_singler_annotations_file "${singler_dir}/singler_annotations.tsv" \
        --output_singler_results_file "${singler_dir}/singler_results.rds" \
        --seed ${params.seed} \
        --threads ${task.cpus}

      # write out meta file
      echo "${Utils.makeJson(meta)}" > "${singler_dir}/scpca-meta.json"
      """
    stub:
      singler_dir = file(meta.singler_dir).name
      """
      mkdir "${singler_dir}"
      echo "${Utils.makeJson(meta)}" > "${singler_dir}/scpca-meta.json"
      """
}


process classify_cellassign {
  container params.SCPCATOOLS_CONTAINER
    publishDir (
      path: "${meta.celltype_publish_dir}",
      mode: 'copy',
      pattern: "${cellassign_dir}"
    )
  label 'mem_32'
  label 'cpus_12'
  input:
    tuple val(meta), path(processed_rds), path(cellassign_reference_file)
  output:
    tuple val(meta.library_id), path(cellassign_dir)
  script:
    cellassign_dir = file(meta.cellassign_dir).name

    """
    # create output directory
    mkdir "${cellassign_dir}"
    
    # Convert SCE to AnnData
    sce_to_anndata.R \
        --input_sce_file "${processed_rds}" \
        --output_rna_h5 processed.hdf5 
        
    # Run CellAssign
    predict_cellassign.py \
      --input_hdf5_file processed.hdf5 
      --output_predictions "${cellassign_dir}/cellassign_predictions.tsv" \
      --reference "${cellassign_reference_file}" \
      --seed ${params.seed} \
      --threads ${task.cpus}
    
    # write out meta file
    echo "${Utils.makeJson(meta)}" > "${cellassign_dir}/scpca-meta.json"
    """
  stub:
    cellassign_dir = file(meta.cellassign_dir).name
    """
    mkdir "${cellassign_dir}"
    echo "${Utils.makeJson(meta)}" > "${cellassign_dir}/scpca-meta.json"
    """
}

// TODO: overhaul this process next
process add_celltypes_to_sce {
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

empty_file = "${projectDir}/assets/NO_FILE"

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
         Utils.parseNA(it.singler_ref_file) ? "${params.singler_models_dir}/${it.singler_ref_file}" : null,
         // cellassign reference file
         Utils.parseNA(it.cellassign_ref_file) ? "${params.cellassign_ref_dir}/${it.cellassign_ref_file}" : null
        ]}

      // create input for typing: [augmented meta, processed_sce]
      celltype_input_ch = processed_sce_channel
        .map{[it[0]["project_id"]] + it}
        .combine(celltype_ch, by: 0)
        // current contents: [project_id, meta, processed_sce, singler_model_file, cellassign_reference_file]
        // add values to meta for later use
        .map{ project_id, meta, processed_sce, singler_model_file, cellassign_reference_file ->
          meta.celltype_publish_dir = "${params.checkpoints_dir}/celltype/${meta.library_id}";
          meta.singler_dir = "${meta.celltype_publish_dir}/${meta.library_id}_singler";
          meta.cellassign_dir = "${meta.celltype_publish_dir}/${meta.library_id}_cellassign";
          meta.singler_model_file = singler_model_file;
          meta.cellassign_reference_file = cellassign_reference_file;
          // return simplified input:
          [meta, processed_sce]
        }

      
      // creates [meta, processed sce, singler model file]
      singler_input_ch = celltype_input_ch
        // add in singler model or empty file
        .map{it.toList() + [file(it[0].singler_model_file ?: empty_file)]}
        // skip if no singleR model file
        .branch{
          missing_ref: it[2].name == "NO_FILE"
          do_singler: true
        }
      

      // perform singleR celltyping and export results
      classify_singler(singler_input_ch.do_singler)
      // singleR output channel: [library_id, singler_results]
      singler_output_ch = singler_input_ch.missing_ref
        .map{[it[0]["library_id"], file(empty_file)]}
        // add in channel outputs
        .mix(classify_singler.out)
      
      // create cellassign input channel: [meta, processed sce, cellassign reference file]
       cellassign_input_ch = celltype_input_ch
        // add in cellassign reference
        .map{it.toList() + [file(it[0].cellassign_reference_file ?: empty_file)]}
        // skip if no cellassign reference file or reference name is not defined
        .branch{
          missing_ref: it[2].name == "NO_FILE"
          do_cellassign: true
        }     

  
      // perform CellAssign celltyping and export results
      classify_cellassign(cellassign_input_ch.do_cellassign)
  
      // cellassign output channel: [library_id, cellassign_dir]
      cellassign_output_ch = cellassign_input_ch.missing_ref
        .map{[it[0]["library_id"], file(empty_file)]}
        // add in channel outputs
        .mix(classify_cellassign.out) 
      
      // prepare input for process to add celltypes to the processed SCE
      assignment_input_ch = processed_sce_channel
        .map{[it[0]["library_id"]] + it}
        // add in singler results
        .join(singler_output_ch, by: 0, failOnMismatch: true, failOnDuplicate: true)
        // add in cell assign results
        .join(cellassign_output_ch, by: 0, failOnMismatch: true, failOnDuplicate: true)
        .map{it.drop(1)} // remove library_id


      // Next PR:
      //add_celltypes_to_sce(assignment_input_ch)
    
      // add back in the unchanged sce files
      // TODO update below with output channel results:
      // export_channel = processed_sce_channel
      //   .map{[it[0]["library_id"]] + it}
      //   // add in unfiltered and filtered sce files
      //   .join(sce_files_channel.map{[it[0]["library_id"], it[1], it[2]]},
      //         by: 0, failOnMismatch: true, failOnDuplicate: true)
      //   // rearrange to be [meta, unfiltered, filtered, processed]
      //   .map{library_id, meta, processed_sce, unfiltered_sce, filtered_sce ->
      //       [meta, unfiltered_sce, filtered_sce, processed_sce]}

    emit: sce_files_channel

}
