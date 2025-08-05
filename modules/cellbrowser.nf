process cellbrowser_site {
  container "${params.CELLBROWSER_CONTAINER}"
  publishDir "${params.results_dir}", mode: 'copy'
  input:
    path(project_metadata)
  output:
    path("scpca_cellbrowser")
  script:
    """
    mkdir -p scpca_cellbrowser
    """
  stub:
    """
    mkdir -p scpca_cellbrowser
    """
}

process cellbrowser_project {
  container "${params.CELLBROWSER_CONTAINER}"
  publishDir "${params.results_dir}", mode: 'copy'
  tag "${project_id}"
  input:
    val(project_id)
    path(project_metadata)
    path(site_dir)
  output:
    tuple val(project_id), path("${site_dir}/${project_id}")
  script:
    """
    mkdir -p "${site_dir}/${project_id}"
    """
  stub:
    """
    mkdir -p "${site_dir}/${project_id}"
    touch "${site_dir}/${project_id}/project_file.html"
    """

}


process cellbrowser_library {
  container "${params.CELLBROWSER_CONTAINER}"
  publishDir "${params.results_dir}/scpca_cellbrowser", mode: 'copy'
  tag "${meta.library_id}"

  input:
    tuple val(meta), path(h5ad_file), path(project_dir)

  output:
    tuple val(meta), path("${project_dir}/${meta.library_id}")

  script:
    """
    mkdir -p "${project_dir}/${meta.library_id}"
    """
  stub:
    """
    mkdir -p "${project_dir}/${meta.library_id}"
    touch "${project_dir}/${meta.library_id}/library_file.html"
    """
}


workflow make_cellbrowser {
  take:
    processed_anndata_ch // channel of tuples [meta, processed_h5ad_file(s)]
  main:
    project_metadata = file(params.project_metafile)
    // export processed anndata files for cellbrowser
    site_base = cellbrowser_site(project_metadata)

    // make a channel of project ids
    project_ch = processed_anndata_ch
      .map{meta, _h5ad -> meta.project_id}
      .unique()

    cellbrowser_project(project_ch, project_metadata, site_base)

    library_ch = processed_anndata_ch
      .map{meta, h5ad_file -> [meta.project_id, meta, h5ad_file]}
      .join(cellbrowser_project.out, by: 0)
      .map{project_id, meta, h5ad_file, project_dir -> [meta, h5ad_file, project_dir]}

    cellbrowser_library(library_ch)

  emit:
    cellbrowser_library.out

}
