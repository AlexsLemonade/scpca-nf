
// process for converting rds files containing an SCE to h5 containing anndata containing the RNA data
process export_anndata{
    container params.SCPCATOOLS_CONTAINER
    label 'mem_16'
    tag "${meta.library_id}"
    publishDir "${params.results_dir}/${meta.project_id}/${meta.sample_id}", mode: 'copy'
    input:
      tuple val(meta), path(sce_file), val(file_type)
    output:
      tuple val(meta), path(rna_hdf5_file), emit: rna_anndata
      tuple val(meta), path(feature_hdf5_file), emit: feature_anndata, optional:true
    script:
      rna_hdf5_file = "${meta.library_id}_${file_type}.hdf5"
      feature_hdf5_file = "${meta.library_id}_${file_type}_feature.hdf5"
      feature_present = meta.feature_type in ["adt", "cellhash"]
      """
      sce_to_anndata.R \
        --input_sce_file ${sce_file} \
        --output_rna_h5 ${rna_hdf5_file} \
        --output_feat_h5 ${feature_hdf5_file} \
        ${feature_present ? "--feature_name ${meta.feature_type}" : ''}
      """
    stub:
      rna_hdf5_file = "${meta.library_id}_${file_type}.hdf5"
      feature_hdf5_file = "${meta.library_id}_${file_type}_feature.hdf5"
      """
      touch ${rna_hdf5_file}
      touch ${feature_hdf5_file}
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

      // combine all anndata files by library id
      // creates anndata channel with [library_id, unfiltered, filtered, processed]
      anndata_ch = export_anndata.out.rna_anndata.mix(export_anndata.out.feature_anndata)
        .map{ meta, hdf5_file -> tuple(
          groupKey(meta.library_id, meta.feature_type in ["adt", "cellhash"]? 6 : 3),
          meta, 
          hdf5_file)
        }
        .groupTuple(by: 0, size: 6, remainder: true)
        .map{ [it[1][0]] +  it[2] }

    emit: anndata_ch

}
