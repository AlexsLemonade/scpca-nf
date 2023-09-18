
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
  label 'mem_8'
  tag "${meta.library_id}"
  publishDir "${params.results_dir}/${meta.project_id}/${meta.sample_id}", mode: 'copy'
  input:
    tuple val(meta), path(processed_hdf5), val(file_type)
  output:
    tuple val(meta), path(processed_hdf5), val(file_type)
  script:
    """
    for file in ${processed_hdf5}; do
      move_counts_anndata.py \
        --anndata_file \${file}
    done
    """
    stub:
       """
       # nothing to do since files don't move
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

     anndata_ch = export_anndata.out
        .branch{
          processed: it[2] == "processed"
          other: true
        }

      // move any normalized counts to X in AnnData
      move_normalized_counts(anndata_ch.processed)

      // PENDING ACTION BELOW:
      // combine all anndata files by library id
      anndata_ch = anndata_ch.other.mix(move_normalized_counts.out)
        // mix with output from moving counts
        .mix(move_normalized_counts.out)
        .map{ meta, hdf5_files, file_type -> tuple(
          meta.library_id, // pull out library id for grouping
          meta,
          hdf5_files // either rna.hdf5 or [ rna.hdf5, feature.hdf5 ]
        )}
        // group by library id result is
        // [library id, [meta, meta, meta], [hdf5 files]]
        .groupTuple(by: 0, size: 3, remainder: true)
        // pull out just 1 meta object and hdf5 files
        // [meta, [hdf5 files]]
        .map{ [it[1][0]] +  it[2] }

    emit: move_normalized_counts.out

}
