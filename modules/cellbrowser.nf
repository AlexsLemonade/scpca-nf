
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
    cellbrowser_config.py \
      --conf_type library\
      --ids "${meta.library_id}"
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
    path site_conf_dir, name: 'site_conf_dir'
  output:
    path "cellbrowser"
  script:
    """
    mkdir -p cellbrowser
    cp -f site_conf_dir/* cellbrowser/

    # create the project config files
    cellbrowser_config.py \
      --outdir cellbrowser \
      --conf_type project \
      --ids ${project_ids.unique(false).join(",")} \
      --project-metadata ${project_metadata}

    # move library directories into place

    for i in \${!libraries[@]}; do
      mkdir -p "cellbrowser/\${projects[\$i]}"
      mv \${libraries[\$i]} "cellbrowser/\${projects[\$i]}/\${libraries[\$i]}"
    done


    """
  stub:
    """
    mkdir -p cellbrowser
    cp -f site_conf_dir/* cellbrowser/

    for project in ${project_ids.unique(false).join(" ")} ; do
      mkdir -p "cellbrowser/\${project}"
      touch "cellbrowser/\${project}/cellbrowser.conf"
      touch "cellbrowser/\${project}/desc.conf"
    done

    libraries=(${library_dirs.join(" ")})
    projects=(${project_ids.join(" ")})
    for i in \${!libraries[@]}; do
      mkdir -p "cellbrowser/\${projects[\$i]}"
      mv \${libraries[\$i]} "cellbrowser/\${projects[\$i]}/\${libraries[\$i]}"
    done
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
