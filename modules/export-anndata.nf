
// process for converting rds files containing sce to h5 containing anndata
process export_anndata{
    container params.SCPCATOOLS_CONTAINER
    label 'mem_16'
    tag "${meta.library_id}"
    publishDir "${params.results_dir}/${meta.project_id}/${meta.sample_id}", mode: 'copy'
    input:
      tuple val(meta), path(sce_file), val(file_type)
    output:
      tuple val(meta), path(hdf5_file)
    script:
      hdf5_file = "${meta.library_id}_${file_type}.hdf5"

      """
      sce_to_anndata.R \
        --input_sce_file ${sce_file} \
        --output_h5_file ${hdf5_file}
      """
    stub:
      hdf5_file = "${meta.library_id}_${file_type}.hdf5"
      """
      touch ${hdf5_file}
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
      // creates anndata channel with [meta, unfiltered, filtered, processed]
      anndata_ch = export_anndata.out
        .map{ [it[0]["library_id"], it[0], it[1]] }
        .groupTuple(by: 0, size: 3, remainder: true)
        .map{ [it[1][0]] +  it[2] }

    emit: anndata_ch

}
