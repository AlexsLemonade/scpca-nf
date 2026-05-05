
include { pullthroughContainer } from '../lib/utils.nf'

process cellbrowser_library {
  container "${pullthroughContainer(params.cellbrowser_container, params.pullthrough_registry)}"
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

    # Rename columns so cbImportScanpy uses the gene_symbols branch in geneStringsFromVar:
    # with gene_ids absent and gene_symbols present it writes features.tsv.gz as ENSG<tab>symbol.
    python3 << 'EOF'
import anndata
adata = anndata.read_h5ad("${h5ad_file}")
adata.var = adata.var.rename(columns={"gene_ids": "ensembl_id"})
if "gene_symbol" in adata.var.columns:
    adata.var["gene_symbols"] = [
        symbol if isinstance(symbol, str) else ensembl_id
        for symbol, ensembl_id in zip(adata.var["gene_symbol"], adata.var["ensembl_id"])
    ]
adata.uns = {}  # remove uns to prevent export error
adata.write_h5ad("cb_input.h5ad")
EOF

    # import data to library directory
    # --proc flag uses the processed expression matrix (lognormalized)
    cbImportScanpy -i "cb_input.h5ad" -o "${meta.library_id}" --clusterField="${params.cellbrowser_default_label}" --proc

    # remove the h5ad files from export
    rm -f "${meta.library_id}"/*.h5ad

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
  container "${pullthroughContainer(params.cellbrowser_container, params.pullthrough_registry)}"
  publishDir "${params.outdir}"
  stageInMode 'copy'
  input:
    tuple val(project_ids), path(library_dirs)
    path project_metadata
    path template_dir, stageAs: "cb_data"
  output:
    path params.cellbrowser_dirname
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
      # move files into place for cellbrowser build
      library_id=\$(basename \${library_dirs[\$i]})
      project_dir="cb_data/\${project_ids[\$i]}/"
      mkdir -p "\${project_dir}"
      mv \${library_dirs[\$i]} "\${project_dir}/\${library_id}"

      # annotate the marker file and if it succeeds, replace the unannotated version with the annotated version.
      cbMarkerAnnotate "\${project_dir}/\${library_id}/markers.tsv" "\${project_dir}/\${library_id}/markers_annotated.tsv"
      if [ -f "\${project_dir}/\${library_id}/markers_annotated.tsv" ]; then
        sed 's/markers.tsv/markers_annotated.tsv/g' -i "\${project_dir}/\${library_id}/cellbrowser.conf"
      fi
    done

    # build the site
    CBDATAROOT=cb_data cbBuild -r -i cb_data/cellbrowser.conf \
      --outDir ${params.cellbrowser_dirname}

    # fix the marker files to be complete (annotated gene symbols were lost in the cbBuild process)
    cellbrowser_rewritemarkers.py \
      --library-ids "${library_dirs.collect{d -> d.name}.join(",")}" \
      --project-ids "${project_ids.join(",")}" \
      --outdir ${params.cellbrowser_dirname}
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

    // use existing output directory if it exists (This caused trouble with rebuilds; we may want to revisit it later)
    // def cb_outdir = file("${params.outdir}/${params.cellbrowser_dirname}", type: 'dir')
    // if (!cb_outdir.exists()) {
    //   cb_outdir.mkdirs()
    // }

    // create single channel of [[project_ids], [library_dirs]]
    project_libs_ch = cellbrowser_library.out
    // only include libraries with umap
     .filter{ it -> it[2] == "true" }
     // use dummy value to group everything together into tuples
     .map{ meta, library_dir, _has_umap ->
        [1, meta.project_id, library_dir]
      }
     .groupTuple()
     .map{ it -> it.drop(1) } // drop the dummy value

    // export processed anndata files for cellbrowser
    cellbrowser_site(
      project_libs_ch,
      file(params.project_metafile),
      file(params.cellbrowser_template_dir, type: 'dir', checkIfExists: true)
    )

  emit:
    cellbrowser_site.out
}
