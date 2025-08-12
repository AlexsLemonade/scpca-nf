process cellbrowser_site {
  container "${params.CELLBROWSER_CONTAINER}"
  publishDir params.cellbrowser_dir, mode: 'copy'
  input:
    val project_ids // list of project ids to create sites for
    path project_metadata
  output:
    path "*"

  script:
    """
    touch cellbrowser.conf
    touch desc.conf
    touch abstract.html


    projects=(${project_ids.join(" ")})
    for project in \${projects[@]}; do
      mkdir -p \${project}
      touch \${project}/project_file.html
    done

    """
  stub:
    """
    touch cellbrowser.conf
    touch desc.conf
    touch abstract.html

    projects=(${project_ids.join(" ")})
    for project in \${projects[@]}; do
      mkdir -p \${project}
      touch \${project}/project_file.html
    done

    """
}


process cellbrowser_library {
  container "${params.CELLBROWSER_CONTAINER}"
  publishDir params.cellbrowser_dir, mode: 'copy'
  tag "${meta.library_id}"

  input:
    tuple val(meta), path(h5ad_file)
    path(site_files)

  output:
    tuple val(meta), path("${meta.project_id}/${meta.library_id}")

  script:
    """
    mkdir -p "${meta.project_id}/${meta.library_id}"
    touch "${meta.project_id}/${meta.library_id}/library_file.html"
    """
  stub:
    """
    mkdir -p "${meta.project_id}/${meta.library_id}"
    touch "${meta.project_id}/${meta.library_id}/library_file.html"
    """
}


workflow build_cellbrowser {
  take:
    processed_anndata_ch // channel of tuples [meta, processed_h5ad_file(s)]
  main:
    project_metadata = file(params.project_metafile)
    // get list of projects
    project_ids = processed_anndata_ch
      .map{it[0].project_id}
      .unique()
      .toList()

    // export processed anndata files for cellbrowser
    cellbrowser_site(project_ids, project_metadata)

    cellbrowser_library(processed_anndata_ch, cellbrowser_site.out)

  emit:
    cellbrowser_library.out

}
