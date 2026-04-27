
include { pullthroughContainer } from '../lib/utils.nf'
// utility process to convert anndata files to h5ad files for downstream processing
// (simpler and fewer options than export_anndata, which is used for exporting files for external use)

process convert_anndata {
  container "${pullthroughContainer(params.scpcatools_anndata_container, params.pullthrough_registry)}"
  label 'mem_16'
  tag "${meta.unique_id}"
  input:
    tuple val(meta), path(processed_rds)
  output:
    tuple val(meta), path("processed.h5ad")
  script:
    """
    sce_to_anndata.R \
        --input_sce_file "${processed_rds}" \
        --output_rna_h5 "processed.h5ad"

    # if the h5ad file was not created successfully, create an empty file so downstream processes don't break
    if [[ ! -e "processed.h5ad" ]]; then
      touch "processed.h5ad"
    fi
    """
  stub:
    """
    # non-empty file for later tests
    echo "processed" > processed.h5ad
    """
}
