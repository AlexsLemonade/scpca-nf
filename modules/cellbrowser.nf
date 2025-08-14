
process cellbrowser_library {
  container "${params.CELLBROWSER_CONTAINER}"
  tag "${meta.library_id}"

  input:
    tuple val(meta), path(h5ad_file)

  output:
    tuple val(meta), path("${meta.library_id}")

  script:
    """
    infile="${meta.library_id}_processed_rna.h5ad" # Name the file in case there are multiple modalities

    # create the library config files
    cellbrowser_config.py \
      --conf_type library \
      --ids "${meta.library_id}" \
      --label-field cluster \
      --sample-ids "${meta.sample_id}"

    # import data to library directory
    cbImportScanpy -i "\${infile}" -o "${meta.library_id}" --clusterField="cluster"

    # remove the h5ad as we won't use it
    rm "${meta.library_id}"/*_processed_rna.h5ad
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
    path "cb_data"
  output:
    path "cb_site"
  script:
    """
    # create the project config files
    cellbrowser_config.py \
      --outdir cb_data \
      --conf_type project \
      --ids ${project_ids.unique(false).join(",")} \
      --project-metadata ${project_metadata}

    # move library directories into place
    library_dirs=(${library_dirs.join(" ")})
    project_ids=(${project_ids.join(" ")})
    for i in \${!library_dirs[@]}; do
      library_id=\$(basename \${library_dirs[\$i]})
      mv \${library_dirs[\$i]} "cb_data/\${project_ids[\$i]}/\${library_id}"
    done

    # build the site
    CBDATAROOT=cb_data cbBuild -r -i cb_data/cellbrowser.conf -o cb_site
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
