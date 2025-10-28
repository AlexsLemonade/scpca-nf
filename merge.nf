#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Workflow to merge SCE objects into a single object.
// This workflow does NOT perform integration, i.e. batch correction.

// parameter check
def check_parameters() {
  def param_error = false

  // check that at least one project has been provided
  if (!params.project) {
    log.error("At least one 'project' must be specified for merging.")
    param_error = true
  }

  // check for provided run file
  if (!file(params.run_metafile).exists()) {
    log.error("The 'run_metafile' file '${params.run_metafile}' can not be found.")
    param_error = true
  }

  if (param_error) {
    System.exit(1)
  }
}


// merge individual SCE objects into one SCE object
process merge_sce {
  container Utils.pullthroughContainer(params.scpcatools_slim_container, params.pullthrough_registry)
  tag "${merge_group_id}"
  label 'mem_max'
  label 'long_running'
  publishDir "${params.results_dir}/${merge_group_id}/merged"
  input:
    tuple val(merge_group_id), val(has_adt), val(library_ids), path(scpca_nf_file)
  output:
    tuple val(merge_group_id), path(merged_sce_file), val(has_adt)
  script:
    input_library_ids = library_ids.join(',')
    input_sces = scpca_nf_file.join(',')
    merged_sce_file = "${merge_group_id}_merged.rds"
    """
    merge_sces.R \
      --input_library_ids "${input_library_ids}" \
      --input_sce_files "${input_sces}" \
      --output_sce_file "${merged_sce_file}" \
      --n_hvg ${params.num_hvg} \
      ${has_adt ? "--include_altexp" : ''} \
      --threads ${task.cpus}
    """
  stub:
    merged_sce_file = "${merge_group_id}_merged.rds"
    """
    touch ${merged_sce_file}
    """

}

// create merge report
process generate_merge_report {
  container Utils.pullthroughContainer(params.scpcatools_reports_container, params.pullthrough_registry)
  tag "${merge_group_id}"
  publishDir "${params.results_dir}/${merge_group_id}/merged"
  label 'mem_max'
  input:
    tuple val(merge_group_id), path(merged_sce_file), val(has_adt)
    path(report_template)
  output:
    path(merge_report)
  script:
    merge_report = "${merge_group_id}_merged-summary-report.html"
    """
    Rscript -e "rmarkdown::render( \
      '${report_template}', \
      output_file = '${merge_report}', \
      params = list(merge_group = '${merge_group_id}', \
                    merged_sce_file = '${merged_sce_file}', \
                    batch_column = 'library_id') \
      )"
    """
  stub:
    merge_report = "${merge_group_id}_merged-summary-report.html"
    """
    touch ${merge_report}
    """
}

process export_anndata {
  container Utils.pullthroughContainer(params.scpcatools_anndata_container, params.pullthrough_registry)
  label 'mem_max'
  label 'long_running'
  tag "${merge_group_id}"
  publishDir "${params.results_dir}/${merge_group_id}/merged", mode: 'copy'
  input:
    tuple val(merge_group_id), path(merged_sce_file), val(has_adt)
  output:
    tuple val(merge_group_id), path("${merge_group_id}_merged_*.h5ad")
  script:
    rna_h5ad_file = "${merge_group_id}_merged_rna.h5ad"
    feature_h5ad_file = "${merge_group_id}_merged_adt.h5ad"
    """
    sce_to_anndata.R \
      --input_sce_file ${merged_sce_file} \
      --output_rna_h5 ${rna_h5ad_file} \
      --output_feature_h5 ${feature_h5ad_file} \
      --is_merged \
      ${has_adt ? "--feature_name adt" : ''}

    # move normalized counts to X in AnnData
    reformat_anndata.py --anndata_file ${rna_h5ad_file} --hvg_name "merged_highly_variable_genes"
    ${has_adt ? "reformat_anndata.py --anndata_file ${feature_h5ad_file} --hvg_name 'none' " : ''}
    """
  stub:
    rna_h5ad_file = "${merge_group_id}_merged_rna.h5ad"
    feature_h5ad_file = "${merge_group_id}_merged_adt.h5ad"
    """
    touch ${rna_h5ad_file}
    ${has_adt ? "touch ${feature_h5ad_file}" : ''}
    """
}

workflow {
  check_parameters()

  // define path to merge template
  def merge_template = "${projectDir}/templates/merge-report.rmd"


  // grab project ids to run
  def project_ids = params.project?.tokenize(',') ?: []

  // grab run ids to include
  def run_ids = params.merge_run_ids?.tokenize(',') ?: []
  // if no run ids, run all
  def run_all = run_ids[0] == "All"

  // read in run metafile and filter to projects of interest
  libraries_ch = channel.fromPath(params.run_metafile)
    .splitCsv(header: true, sep: '\t')
    // filter to only include specified project ids
    .filter{ it.scpca_project_id in project_ids }
    // filter to run all ids or just specified ones
    .filter{ it ->
      run_all
      || (it.scpca_run_id in run_ids)
      || (it.scpca_library_id in run_ids)
      || (it.scpca_sample_id in run_ids)
    }
    .map{ it ->
      [
        project_id: it.scpca_project_id,
        library_id: it.scpca_library_id,
        sample_id: it.scpca_sample_id.split(";").sort().join(","),
        seq_unit: it.seq_unit,
        technology: it.technology
      ]
    }

  // get all projects that contain at least one library with CITEseq
  adt_projects = libraries_ch
    .filter{ it.technology.startsWith('citeseq') }
    .collect{ it.project_id }
    .map{ it -> it.unique() }

  multiplex_projects = libraries_ch
    .filter{ it.technology.startsWith('cellhash') }
    .collect{ it.project_id }
    .map{ it -> it.unique() }

  oversized_projects = libraries_ch
    .filter{ it.technology.startsWith("10x") } // only count single-cell or single-nuclei libraries, no cell hash, ADT, bulk or spatial
    // pull out project id for grouping
    .map{ it -> [it.project_id, it] }
    .groupTuple(by: 0) // group by project id
    .filter{ it[1].size() > params.max_merge_libraries } // get projects with more samples than max merge
    .collect{ it[0] } // in this specific case, it[0] is the project_id

  filtered_libraries_ch = libraries_ch
    // only include single-cell/single-nuclei which ensures we don't try to merge libraries from spatial or bulk data
    .filter{ it.seq_unit in ['cell', 'nucleus'] }
    // remove any multiplexed projects or oversized projects
    // future TODO: only filter library ids that are multiplexed, but keep all other non-multiplexed libraries
    .branch{ it ->
      multiplexed: it.project_id in multiplex_projects.getVal()
      oversized: it.project_id in oversized_projects.getVal()
      single_sample: true
    }

  filtered_libraries_ch.multiplexed
    .unique{ it.project_id }
    .subscribe{ it ->
      log.warn("Not merging ${it.project_id} because it contains multiplexed libraries.")
    }

  filtered_libraries_ch.oversized
    .unique{ it.project_id }
    .subscribe{ it ->
      log.warn("Not merging ${it.project_id} because it contains too many libraries.")
    }

  // print out warning message for any libraries not included in merging
  filtered_libraries_ch.single_sample
    .map{ it ->
      def processed = file("${params.results_dir}/${it.project_id}/${it.sample_id}/${it.library_id}_processed.rds")
      def meta_json = file("${params.results_dir}/${it.project_id}/${it.sample_id}/${it.library_id}_metadata.json")
      [it.library_id, processed, meta_json ]
    }
    .subscribe{ library_id, processed, meta_json ->
      if (!processed.exists() || !(processed.size() > 0 || library_id.startsWith("STUBL"))) {
        log.warn("Processed files do not exist for ${library_id}. This library will not be included in the merged object.")
      }
      else if (!(meta_json.exists() && meta_json.size() > 0)) {
        log.warn("Metadata file does not exist for ${library_id}. This library will not be included in the merged object.")
      }
      else if (Utils.getMetaVal(meta_json, "processed_cells") < 3) {
        log.warn("Library ${library_id} has fewer than 3 cells. This library will not be included in the merged object.")
      }
    }

  grouped_libraries_ch = filtered_libraries_ch.single_sample
    .map{ it ->
      def processed = file("${params.results_dir}/${it.project_id}/${it.sample_id}/${it.library_id}_processed.rds")
      def meta_json = file("${params.results_dir}/${it.project_id}/${it.sample_id}/${it.library_id}_metadata.json")
      [it.project_id, it.library_id, processed, meta_json]
    }
    // only include libraries that have been processed through scpca-nf and have at least 3 cells
    .filter{ _project_id, library_id, processed, meta_json ->
      (processed.exists() && processed.size() > 0 && Utils.getMetaVal(meta_json, "processed_cells") >= 3)
      || library_id.startsWith("STUBL")
    }
    // remove metadata file
    .map{ project_id, library_id, processed, _meta_json ->
      [project_id, library_id, processed]
    }
    // only one row per library ID, this removes all the duplicates that may be present due to CITE/hashing
    .unique()
    // group tuple by project id: [project_id, [library_id1, library_id2, ...], [sce_file1, sce_file2, ...]]
    .groupTuple(by: 0)
    // add in boolean for if project contains samples with adt
    .map{ project_id, library_id_list, sce_file_list ->
      [project_id, project_id in adt_projects.getVal(), library_id_list, sce_file_list]
    }
    .branch{ it ->
      def merged_sce = file("${params.results_dir}/${it[0]}/merged/${it[0]}_merged.rds")
      has_merge: merged_sce.exists() && params.reuse_merge
      make_merge: true
    }

  pre_merged_ch = grouped_libraries_ch.has_merge
    .map{ project_id, has_adt, _library_id_list, _sce_file_list ->
      def merged_file = file("${params.results_dir}/${project_id}/merged/${project_id}_merged.rds")
      [project_id, merged_file, has_adt]
    }

  // merge SCE objects
  merge_sce(grouped_libraries_ch.make_merge)

  merged_ch = merge_sce.out.mix(pre_merged_ch)


  // generate merge report
  generate_merge_report(merged_ch, file(merge_template))

  // export merged objects to AnnData
  export_anndata(merged_ch)
}
