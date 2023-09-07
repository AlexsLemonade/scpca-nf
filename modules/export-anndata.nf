
// process for converting rds files containing an SCE to h5 containing anndata containing the RNA data
process export_anndata{
    container params.SCPCATOOLS_CONTAINER
    label 'mem_16'
    tag "${meta.library_id}"
    publishDir "${params.results_dir}/${meta.project_id}/${meta.sample_id}", mode: 'copy'
    input:
      tuple val(meta), path(sce_file), val(file_type)
    output:
      tuple val(meta), path("${meta.library_id}_${file_type}*.hdf5"), val(file_type)
    script:
      rna_hdf5_file = "${meta.library_id}_${file_type}_rna.hdf5"
      feature_hdf5_file = "${meta.library_id}_${file_type}_${meta.feature_type}.hdf5"
      feature_present = meta.feature_type in ["adt", "cellhash"]
      """
      sce_to_anndata.R \
        --input_sce_file ${sce_file} \
        --output_rna_h5 ${rna_hdf5_file} \
        --output_feature_h5 ${feature_hdf5_file} \
        ${feature_present ? "--feature_name ${meta.feature_type}" : ''}
      """
    stub:
      rna_hdf5_file = "${meta.library_id}_${file_type}_rna.hdf5"
      feature_hdf5_file = "${meta.library_id}_${file_type}_${meta.feature_type}.hdf5"
      """
      touch ${rna_hdf5_file}
      touch ${feature_hdf5_file}
      """
}

process move_normalized_counts{
  container params.SCPCATOOLS_CONTAINER
  label 'mem_4'
  tag "${meta.library_id}"
  publishDir "${params.results_dir}/${meta.project_id}/${meta.sample_id}", mode: 'copy'
  input:
    tuple val(meta), path(processed_hdf5), val(file_type)
  output:
    tuple val(meta), path(processed_hdf5), val(file_type)
  script:
    input_hdf5_files = processed_hdf5.join(',')
    """
    move_counts_anndata.py \
      --input_hdf5_file ${input_hdf5_files}
    """
}


workflow sce_to_anndata{
    take:
      sce_files_ch
    main:
      sce_ch = sce_files_ch
        // make tuple of [meta, sce_file, type of file]
        .flatMap{[[it[0], it[1], "unfiltered"],
                  [it[0], it[2], "filtered"],
                  [it[0], it[3], "processed"]
                 ]}

      // export each anndata file
      export_anndata(sce_ch)

      processed_anndata_ch = export_anndata.out
        .filter{ it[2] == "processed"}

      // move any normalized counts to X in AnnData
      move_normalized_counts(processed_anndata_ch)

      // combine all anndata files by library id
      // creates anndata channel with [library_id, unfiltered, filtered, processed]
      anndata_ch = export_anndata.out
        .join(move_normalized_counts.out, by: [0-2])
        .map{ meta, hdf5_files, file_type -> tuple(
          meta.library_id,
          meta,
          hdf5_files
        )}
        .groupTuple(by: 0, size: 3, remainder: true)
        .map{ [it[1][0]] +  it[2] }

    emit: anndata_ch

}
