
include { pullthroughContainer } from '../lib/utils.nf'

// Convert one processed AnnData file into a Cell Browser input directory.
// Also annotates the marker table
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

    # annotate the marker file and if it succeeds, replace the unannotated version
    # with the annotated version
    cbMarkerAnnotate "${meta.library_id}/markers.tsv" "${meta.library_id}/markers_annotated.tsv"
    if [ -f "${meta.library_id}/markers_annotated.tsv" ]; then
      sed 's/markers.tsv/markers_annotated.tsv/g' -i "${meta.library_id}/cellbrowser.conf"
    fi

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

// Build the config-only skeleton of the dataset hierarchy
process cellbrowser_project {
  container "${pullthroughContainer(params.cellbrowser_container, params.pullthrough_registry)}"

  input:
    val project_ids
    path project_metadata
    path template_dir, stageAs: "cb_template"

  output:
    path "cb_data"

  script:
    """
    # -L because the staged template is a symlink and cellbrowser_config.py
    # writes into this directory
    cp -rL cb_template cb_data

    cellbrowser_config.py \
      --outdir cb_data \
      --conf_type project \
      --ids ${project_ids.unique(false).join(",")} \
      --project-metadata ${project_metadata}
    """
  stub:
    """
    cp -rL cb_template cb_data
    for project_id in ${project_ids.unique(false).join(" ")}; do
      mkdir -p "cb_data/\${project_id}"
      touch "cb_data/\${project_id}/cellbrowser.conf"
      touch "cb_data/\${project_id}/desc.conf"
    done
    """
}

// Convert one library dataset. This is the expensive step
process cellbrowser_dataset {
  container "${pullthroughContainer(params.cellbrowser_container, params.pullthrough_registry)}"
  tag "${meta.library_id}"
  label 'long_running'

  input:
    tuple val(meta), path(library_dir)
    path skeleton, stageAs: "cb_skeleton"

  output:
    tuple val(meta),
          // only output the subdirectory for this library, not the whole site
          path("built/${meta.project_id}/${meta.library_id}"),
          path("cb_input_${meta.library_id}")

  script:
    """
    # assemble a minimal dataRoot: skeleton confs plus this one library
    cp -rL cb_skeleton cb_data
    mkdir -p "cb_data/${meta.project_id}"
    mv ${library_dir} "cb_data/${meta.project_id}/${meta.library_id}"

    # no -r: -r ignores -i entirely and globs **/cellbrowser.conf from the CWD
    CBDATAROOT=cb_data cbBuild \
      -i "cb_data/${meta.project_id}/${meta.library_id}/cellbrowser.conf" \
      -o built

    # cellbrowser_rewritemarkers.py runs in cellbrowser_site and reads
    # cb_data/<project_id>/<library_id>/{markers_annotated.tsv,features.tsv.gz}.
    # Forward just those two files so the site process can reconstruct that
    # layout without staging the whole library input directory.
    library_in="cb_data/${meta.project_id}/${meta.library_id}"
    # keep the input directories unique for later staging
    marker_in="cb_input_${meta.library_id}/${meta.project_id}/${meta.library_id}"
    mkdir -p "\${marker_in}"
    for marker_file in features.tsv.gz markers_annotated.tsv; do
      if [ -f "\${library_in}/\${marker_file}" ]; then
        cp "\${library_in}/\${marker_file}" "\${marker_in}/"
      fi
    done
    """
  stub:
    """
    mkdir -p "built/${meta.project_id}/${meta.library_id}/markers/markers_0"
    touch "built/${meta.project_id}/${meta.library_id}/dataset.json"
    touch "built/${meta.project_id}/${meta.library_id}/cellbrowser.conf"
    mkdir -p "cb_input_${meta.library_id}/${meta.project_id}/${meta.library_id}"
    """
}

// Assemble the site.
process cellbrowser_site {
  container "${pullthroughContainer(params.cellbrowser_container, params.pullthrough_registry)}"
  publishDir "${params.outdir}"

  input:
    tuple val(project_ids), path(dataset_dirs), path(marker_inputs)
    path skeleton, stageAs: "cb_skeleton"

  output:
    path params.cellbrowser_dirname

  script:
    """
    cp -rL cb_skeleton cb_data

    # rebuild the cb_data/<project_id>/<library_id>/ layout that
    # cellbrowser_rewritemarkers.py expects, from the forwarded marker inputs.
    # cp -r merges into the existing project directories from the skeleton.
    for marker_dir in ${marker_inputs.join(" ")}; do
      [ -d "\${marker_dir}" ] || continue
      cp -rL "\${marker_dir}/." cb_data/
    done

    # place the pre-built dataset directories into the output hierarchy.
    mkdir -p ${params.cellbrowser_dirname}
    dataset_dirs=(${dataset_dirs.join(" ")})
    project_ids=(${project_ids.join(" ")})
    for i in \${!dataset_dirs[@]}; do
      library_id=\$(basename \${dataset_dirs[\$i]})
      project_dir="${params.cellbrowser_dirname}/\${project_ids[\$i]}"
      mkdir -p "\${project_dir}"
      # use -L to make sure we are getting a self-contained directory
      cp -rL "\${dataset_dirs[\$i]}" "\${project_dir}/\${library_id}"
    done

    # single index pass: top-level collection plus every project collection.
    # -i is action="append", and non-recursive build() accumulates all of them
    # into one rebuildCollections call.
    conf_args=(-i cb_data/cellbrowser.conf)
    for project_id in \$(printf '%s\\n' "\${project_ids[@]}" | sort -u); do
      conf_args+=(-i "cb_data/\${project_id}/cellbrowser.conf")
    done
    CBDATAROOT=cb_data cbBuild "\${conf_args[@]}" -o ${params.cellbrowser_dirname}

    # fix the marker files to be complete (annotated gene symbols were lost in the cbBuild process)
    cellbrowser_rewritemarkers.py \
      --library-ids "${dataset_dirs.collect{d -> d.name}.join(",")}" \
      --project-ids "${project_ids.join(",")}" \
      --cb-data cb_data \
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

    // only include libraries with umap
    libraries_ch = cellbrowser_library.out
      .filter{ it -> it[2] == "true" }
      .map{ meta, library_dir, _has_umap -> [meta, library_dir] }

    // one config skeleton covering every project present in the run
    project_ids_ch = libraries_ch
      .map{ meta, _library_dir -> meta.project_id }
      .unique()
      .collect(sort: true)

    cellbrowser_project(
      project_ids_ch,
      file(params.project_metafile),
      file(params.cellbrowser_template_dir, type: 'dir', checkIfExists: true)
    )
    // .first() converts the single-element queue channel into a value channel
    // so it can be consumed by every cellbrowser_dataset task
    skeleton_ch = cellbrowser_project.out.first()

    // convert each library in parallel
    cellbrowser_dataset(libraries_ch, skeleton_ch)

    // create single channel of [[project_ids], [dataset_dirs], [marker_inputs]]
    project_datasets_ch = cellbrowser_dataset.out
      // use dummy value to group everything together into tuples
      .map{ meta, dataset_dir, marker_input ->
        [1, meta.project_id, dataset_dir, marker_input]
      }
      .groupTuple()
      .map{ it -> it.drop(1) } // drop the dummy value

    cellbrowser_site(project_datasets_ch, skeleton_ch)

  emit:
    cellbrowser_site.out
}
