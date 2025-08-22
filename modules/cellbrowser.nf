
process cellbrowser_library {
  container "${params.CELLBROWSER_CONTAINER}"
  tag "${meta.library_id}"

  input:
    tuple val(meta), path(h5ad_file, arity: '1')

  output:
    tuple val(meta), path("${meta.library_id}"), env('has_umap')

  script:
    """
    # create the library config files
    cellbrowser_config.py \
      --conf_type library \
      --ids "${meta.library_id}" \
      --label-field ${params.cellbrowser_default_label} \
      --sample-ids "${meta.sample_id}" \
      --h5ad-file "${h5ad_file}"

    # import data to library directory
    cbImportScanpy -i "${h5ad_file}" -o "${meta.library_id}" --clusterField="${params.cellbrowser_default_label}"

    # remove the h5ad from the imported files as we won't use it
    rm "${meta.library_id}"/*_processed_rna.h5ad

    # Check that the umap coordinates were output
    if [ -f "${meta.library_id}/umap_coords.tsv" ]; then
      has_umap="true"
    else
      has_umap="false"
    fi
    """
  stub:
    """
    mkdir -p "${meta.library_id}"
    # selection of expected files
    touch "${meta.library_id}/cellbrowser.conf"
    touch "${meta.library_id}/desc.conf"
    touch "${meta.library_id}/matrix.mtx.gz"
    has_umap="true"
    """
}

process cellbrowser_site {
  container "${params.CELLBROWSER_CONTAINER}"
  publishDir "${params.outdir}"
  input:
    tuple val(project_ids), path(library_dirs)
    path project_metadata
    path "cb_data"
    path "${params.cellbrowser_dirname}"
  output:
    path "${params.cellbrowser_dirname}"
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
    CBDATAROOT=cb_data cbBuild -r -i cb_data/cellbrowser.conf \
      --redo ${params.cellbrowser_rebuild ? "matrix" : "meta"} \
      --outDir ${params.cellbrowser_dirname}
    """
  stub:
    """
    mkdir -p ${params.cellbrowser_dirname}
    touch ${params.cellbrowser_dirname}/index.html
    """
}


workflow cellbrowser_build {
  take:
    processed_anndata_ch // channel of tuples [meta, processed_h5ad_file]
  main:
    cellbrowser_library(processed_anndata_ch)

    // use existing output directory if it exists
    def cb_outdir = file("${params.outdir}/${params.cellbrowser_dirname}", type: 'dir')
    if (!cb_outdir.exists()) {
      cb_outdir.mkdirs()
    }

    // create single channel of [[project_ids], [library_dirs]]
    project_libs_ch = cellbrowser_library.out
    // only include libraries with umap
     .filter{it.has_umap == "true" }
     // use dummy value to group everything together into tuples
     .map{meta, library_dir, _has_umap -> [1, meta.project_id, library_dir] }
     .groupTuple()
     .map{it -> it.drop(1)}

    // export processed anndata files for cellbrowser
    cellbrowser_site(
      project_libs_ch,
      file(params.project_metafile),
      file(params.cellbrowser_template_dir, type: 'dir', checkIfExists: true),
      cb_outdir
    )

  emit:
    cellbrowser_site.out
}
