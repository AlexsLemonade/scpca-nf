
include { getVersions; makeJson; readMeta; getMetaVal; pullthroughContainer } from '../lib/utils.nf'

// process for converting rds files containing an SCE to h5 containing anndata containing the RNA data
process export_anndata {
    container pullthroughContainer(params.scpcatools_anndata_container, params.pullthrough_registry)
    label 'mem_16'
    tag "${meta.unique_id}"
    publishDir "${params.results_dir}/${meta.project_id}/${meta.sample_id}", mode: 'copy'
    input:
      tuple val(meta), path(sce_file), val(file_type)
    output:
      tuple val(meta), path("${meta.unique_id}_${file_type}_*.h5ad"), val(file_type)
    script:
      rna_h5ad_file = "${meta.unique_id}_${file_type}_rna.h5ad"
      feature_h5ad_file = "${meta.unique_id}_${file_type}_${meta.feature_type}.h5ad"
      pca_meta_file = "${meta.unique_id}_${file_type}_pca.tsv"
      feature_present = meta.feature_type in ["adt"]
      """
      sce_to_anndata.R \
        --input_sce_file ${sce_file} \
        --output_rna_h5 ${rna_h5ad_file} \
        --output_feature_h5 ${feature_h5ad_file} \
        --output_pca_tsv ${pca_meta_file} \
        ${feature_present ? "--feature_name ${meta.feature_type}" : ''} \
        ${file_type != "processed" ? "--compress_output" : ''}

    # move any normalized counts to X in AnnData, convert matrices, and add PCA metadata
    if [ "${file_type}" = "processed" ]; then
      reformat_anndata.py --anndata_file ${rna_h5ad_file} --pca_meta_file ${pca_meta_file}
      # move counts in feature data, if the file exists
      if [ -f "${feature_h5ad_file}" ]; then
        reformat_anndata.py --anndata_file ${feature_h5ad_file} --hvg_name "none"
      fi
    fi

      """
    stub:
      rna_h5ad_file = "${meta.unique_id}_${file_type}_rna.h5ad"
      feature_h5ad_file = "${meta.unique_id_prefix}_${file_type}_${meta.feature_type}.h5ad"
      feature_present = meta.feature_type in ["adt"]
      """
      touch ${rna_h5ad_file}
      ${feature_present ? "touch ${feature_h5ad_file}" : ''}
      """
}


workflow sce_to_anndata {
  take:
    // tuple of [meta, unfiltered rds, filtered rds, processed rds, metadata json]
    sce_files_ch
  main:
    sce_ch = sce_files_ch
      // spread files so only one type of file gets passed through to the process at a time
      // make tuple of [meta, sce_file, type of file, metadata.json]
      .flatMap{ meta, unfiltered, filtered, processed, meta_json ->
        [
          [meta, unfiltered, "unfiltered", meta_json],
          [meta, filtered, "filtered", meta_json],
          [meta, processed, "processed", meta_json]
        ]
      }
      // remove any sce files that don't have enough cells in the sce object
      // number of cells are stored in each metadata.json file
      .filter{ _meta, _sce_file, file_type, meta_json ->
        def cells = getMetaVal(file(meta_json), "${file_type}_cells");
        cells == '' || cells > 2  // if no cell count, keep file for testing, otherwise require at least 3 cells
      }
      // remove metadata.json file from tuple
      .map{ it -> it.dropRight(1) }

    // export each anndata file
    export_anndata(sce_ch)

    // get processed anndata files for cellbrowser
    processed_ch = export_anndata.out
      .filter{ _meta, _h5ad_file, file_type ->
        file_type == "processed"
      }
      .map{ it -> it.dropRight(1) } // drop the file type


    // combine all anndata files by unique id
    anndata_ch = export_anndata.out
      .map{ meta, h5ad_files, _file_type ->
        [
          meta.unique_id, // pull out unique id for grouping
          meta,
          h5ad_files // either *_rna.h5ad or [ *_rna.h5ad, *_feature.h5ad ]
        ]
      }
      // group by unique id result is
      // [unique id, [meta, meta, meta], [h5ad files]]
      .groupTuple(by: 0, size: 3, remainder: true)
      // pull out just 1 meta object and h5ad files
      // [meta, [h5ad files]]
      .map{ _unique_id, meta_list, h5ad_files ->
        [meta_list[0], h5ad_files]
      }

  emit:
    complete = anndata_ch
    processed = processed_ch

}
