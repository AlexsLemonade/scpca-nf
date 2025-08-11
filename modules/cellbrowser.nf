process cellbrowser_site {
  container "${params.CELLBROWSER_CONTAINER}"
  publishDir "${params.cellbrowser_dir}", mode: 'copy'
  input:
    val project_ids // list of project ids to create sites for
    path project_metadata
  output:
    path "scpca_cellbrowser"
  script:
    """
    # create the site directory
    mkdir -p scpca_cellbrowser

    for project in \${project_ids}; do
      mkdir -p scpca_cellbrowser/\${project}
      touch scpca_cellbrowser/\${project}/project_file.html
    done

    """
  stub:
    """
    project_dir=\$(basename ${params.cellbrowser_dir})

    for project in \${projects}; do
      mkdir -p \${project_dir}/\${project}
      touch \${project_dir}/\${project}/project_file.html
    done

    """
}


process cellbrowser_library {
  container "${params.CELLBROWSER_CONTAINER}"
  publishDir "${params.results_dir}", mode: 'copy'
  tag "${meta.library_id}"

  input:
    tuple val(meta), path(h5ad_file)
    path(site_base)

  output:
    tuple val(meta), path("${site_base}/${meta.project_id}/${meta.library_id}")

  script:
    """
    mkdir -p "${site_base}/${meta.project_id}/${meta.library_id}"
    """
  stub:
    """
    mkdir -p "${site_base}/${meta.project_id}/${meta.library_id}"
    touch "${site_base}/${meta.project_id}/${meta.library_id}/library_file.html"
    """
}


workflow build_cellbrowser {
  take:
    processed_anndata_ch // channel of tuples [meta, processed_h5ad_file(s)]
  main:
    project_metadata = file(params.project_metafile)
    // get list of projects
    project_ids = processed_anndata_ch
      .map{it[0].scpca_project_id}
      .unique()
      .toList()

    // export processed anndata files for cellbrowser
    site_base = cellbrowser_site(project_ids, project_metadata)

    cellbrowser_library(processed_anndata_ch, site_base)

  emit:
    cellbrowser_library.out

}
