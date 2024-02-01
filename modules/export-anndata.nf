
// process for converting rds files containing an SCE to h5 containing anndata containing the RNA data
process export_anndata{
    container params.SCPCATOOLS_CONTAINER
    label 'mem_16'
    tag "${meta.library_id}"
    publishDir "${params.results_dir}/${meta.project_id}/${meta.sample_id}", mode: 'copy'
    input:
      tuple val(meta), path(sce_file), val(file_type)
    output:
      tuple val(meta), path("${meta.library_id}_${file_type}_*.hdf5"), val(file_type)
    script:
      rna_hdf5_file = "${meta.library_id}_${file_type}_rna.hdf5"
      feature_hdf5_file = "${meta.library_id}_${file_type}_${meta.feature_type}.hdf5"
      feature_present = meta.feature_type in ["adt"]
      """
      sce_to_anndata.R \
        --input_sce_file ${sce_file} \
        --output_rna_h5 ${rna_hdf5_file} \
        ${feature_present ? "--output_feature_h5 ${feature_hdf5_file}" : ''} \
        ${feature_present ? "--feature_name ${meta.feature_type}" : ''}

      # move any normalized counts to X in AnnData
      if [ "${file_type}" = "processed" ]; then
        move_counts_anndata.py --anndata_file ${rna_hdf5_file}
        # move counts in feature data, if the file exists
        if [ -f "${feature_hdf5_file}" ]; then
          move_counts_anndata.py --anndata_file ${feature_hdf5_file}
        fi
      fi

      """
    stub:
      rna_hdf5_file = "${meta.library_id}_${file_type}_rna.hdf5"
      feature_hdf5_file = "${meta.library_id}_${file_type}_${meta.feature_type}.hdf5"
      feature_present = meta.feature_type in ["adt"]
      """
      touch ${rna_hdf5_file}
      ${feature_present ? "touch ${feature_hdf5_file}" : ''}
      """
}


workflow sce_to_anndata{
    take:
      // tuple of [meta, unfiltered rds, filtered rds, processed rds, metadata json]
      sce_files_ch
    main:
      sce_ch = sce_files_ch
        // spread files so only one type of file gets passed through to the process at a time
        // make tuple of [meta, sce_file, type of file, metadata.json]
        .flatMap{[
          [it[0], it[1], "unfiltered", it[4]],
          [it[0], it[2], "filtered", it[4]],
          [it[0], it[3], "processed", it[4]]
        ]}
        // remove any sce files that don't have enough cells in the sce object
        // number of cells are stored in each metadata.json file
        .filter{
          cells = Utils.getMetaVal(file(it[3]), "${it[2]}_cells");
          cells != '' ? cells > 1 : true // if no cell count, keep file (for testing)
        }
        // remove metadata.json file from tuple
        .map{it.dropRight(1)}

      // export each anndata file
      export_anndata(sce_ch)

      // combine all anndata files by library id
      anndata_ch = export_anndata.out
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

    emit: anndata_ch

}
