
process cellbrowser_library {
  container "${params.CELLBROWSER_CONTAINER}"
  tag "${meta.library_id}"

  input:
    tuple val(meta), path(h5ad_file)

  output:
    tuple val(meta), path("${meta.library_id}")

  script:
    """
    mkdir -p "${meta.library_id}"
    # selection of expected files
    touch "${meta.library_id}/cellbrowser.conf"
    touch "${meta.library_id}/desc.conf"
    touch "${meta.library_id}/matrix.mtx.gz"
    """
  stub:
    """
    mkdir -p "${meta.library_id}"
    # selection of expected files
    touch "${meta.library_id}/cellbrowser.conf"
    touch "${meta.library_id}/desc.conf"
    touch "${meta.library_id}/matrix.mtx.gz"
    """
}

process cellbrowser_site {
  container "${params.CELLBROWSER_CONTAINER}"
  publishDir "${params.cellbrowser_dir}"
  input:
    tuple val(project_ids), path(library_dirs)
    path project_metadata
    path site_conf_dir
  output:
    path "cb_site"
  script:
    """
    mkdir cb_site
    touch cb_site/index.html
    """
  stub:
    """
    mkdir cb_site
    touch cb_site/index.html
    """
}


workflow cellbrowser_build {
  take:
    processed_anndata_ch // channel of tuples [meta, processed_h5ad_file]
  main:
    cellbrowser_library(processed_anndata_ch)

    // create single channel of [[project_ids], [library_dirs]]
    project_libs_ch = cellbrowser_library.out
     // use dummy value to grouping everything together into tuples
     .map{meta, library_dir -> [1, meta.project_id, library_dir] }
     .groupTuple()
     .map{it -> it.drop(1)}

    // export processed anndata files for cellbrowser
    cellbrowser_site(
      project_libs_ch,
      file(params.project_metafile),
      file(params.cellbrowser_template_dir, type: 'dir', checkIfExists: true)
    )

  emit:
    cellbrowser_site.out
}
