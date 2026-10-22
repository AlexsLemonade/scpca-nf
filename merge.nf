#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Include utility functions
include { makeJson; getMetaVal; pullthroughContainer } from './lib/utils.nf'

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
  container "${pullthroughContainer(params.scpcatools_slim_container, params.pullthrough_registry)}"
  tag "${merge_group_id}"
  label 'mem_max'
  label 'long_running'
  publishDir "${params.results_dir}/${merge_group_id}/merged"
  input:
    tuple val(merge_group_id), val(has_adt), val(unique_ids), path(scpca_nf_file)
  output:
    tuple val(merge_group_id), path(merged_sce_file), val(has_adt)
  script:
    // TODO: update the merge script and function to refer to things as unique IDs 
    // that's how they are treated but we use library_ids throughout the code and as the name of the batch_column that gets saved in the merged object
    // I think we want to keep "library_id" since unique_ids don't mean anything on the portal side so changing a column name seems not ideal? 
    input_library_ids = unique_ids.join(',')
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
  container "${pullthroughContainer(params.scpcatools_reports_container, params.pullthrough_registry)}"
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
  container "${pullthroughContainer(params.scpcatools_anndata_container, params.pullthrough_registry)}"
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

process check_sce {
  container "${pullthroughContainer(params.scpcatools_slim_container, params.pullthrough_registry)}"
  label 'mem_16'
  tag "${merge_group_id}"
  input:
    tuple val(merge_group_id), 
          path(merged_sce_file)
    path(reference_sce_file)
  output:
    tuple val(merge_group_id), 
          path("${merge_group_id}_formatting_errors.txt")
  script:
    """
    sce_formatting_checks.R \
        --sce_file ${merged_sce_file} \
        --object_type merged \
        --reference_file ${reference_sce_file} \
        --output_file ${merge_group_id}_formatting_errors.txt
    """
  stub:
    """
    touch ${merge_group_id}_formatting_errors.txt
    """
}

process check_anndata {
  container "${pullthroughContainer(params.scpcatools_anndata_container, params.pullthrough_registry)}"
  label 'mem_16'
  tag "${merge_group_id}"
  input:
    tuple val(merge_group_id), 
          path(anndata_file)
    path(reference_anndata_file)
  output:
    tuple val(merge_group_id), 
          path("${merge_group_id}_formatting_errors.txt")
  script:
    """
    anndata_formatting_checks.py \
        --anndata_file ${anndata_file} \
        --object_type merged \
        --reference_file ${reference_anndata_file} \
        --output_file ${merge_group_id}_formatting_errors.txt
    """
  stub:
    """
    touch ${merge_group_id}_formatting_errors.txt
    """
}

process compile_errors {
  container "${pullthroughContainer(params.scpcatools_slim_container, params.pullthrough_registry)}"
  publishDir "${params.format_checks_dir}", mode: 'copy'
  label 'mem_8'
  input:
   // more than one error file per library (unfiltered, filtered, processed), so use stageAs
   path error_files, stageAs: "error_file_*.txt"
   path output_filename
  output:
    path output_filename
  script:
    """
    # check if all files are empty, if so print a success message
    # otherwise concatenate all error files into one output file
    errors="\$(cat ${error_files})"
    if [ -z "\${errors}" ]; then
      echo "No formatting errors found." > ${output_filename}
    else
      echo "\${errors}" > ${output_filename}
    fi
    """
  stub: 
    """
    touch ${output_filename}
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

  // technologies where each sample in a library should become its own row
  def flex_technologies = ["10xflex_v1.1_multi"]
  // seq units that are included in merging
  def seq_units_for_merge = ['cell', 'nucleus', 'nucleus-FFPE']

  libraries_ch = channel.fromPath(params.run_metafile)
      .splitCsv(header: true, sep: '\t')
      // filter to only include specified project ids
      .filter{ it.scpca_project_id in project_ids }
      .flatMap{ row ->
        def technology = row.technology.toLowerCase()
        def sample_ids = row.scpca_sample_id.split(";").sort()
        def library_meta = [
          project_id: row.scpca_project_id,
          run_id: row.scpca_run_id,
          library_id: row.scpca_library_id,
          seq_unit: row.seq_unit,
          technology: technology
        ]

        if (technology in flex_technologies) {
          // flex: one row per library-sample pair
          return sample_ids.collect{ sample_id ->
            library_meta + [
              sample_id: sample_id,
              unique_id: "${row.scpca_library_id}-${sample_id}".toString()
            ]
          }
        } else {
          // everything else, including cellhash multiplexed libraries: one row per library
          return [
            library_meta + [
              sample_id: sample_ids.join(","),
              unique_id: row.scpca_library_id
            ]
          ]
        }
      }
      // filter to run all ids or just specified ones
      // done after splitting flex libraries so individual sample ids can be matched
      .filter{ it ->
        run_all
        || (it.run_id in run_ids)
        || (it.library_id in run_ids)
        || it.sample_id.split(",").any{ sample_id -> sample_id in run_ids }
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

  filtered_libraries_ch = libraries_ch
    // only include single-cell/single-nuclei which ensures we don't try to merge libraries from spatial or bulk data
    .filter{ it.seq_unit in seq_units_for_merge }
    // remove any multiplexed projects
    // future TODO: only filter library ids that are multiplexed, but keep all other non-multiplexed libraries
    .branch{ it ->
      multiplexed: it.project_id in multiplex_projects.getVal()
      single_sample: true
    }

  filtered_libraries_ch.multiplexed
    .unique{ it.project_id }
    .subscribe{ it ->
      log.warn("Not merging ${it.project_id} because it contains multiplexed libraries.")
    }

  // print out warning message for any libraries not included in merging
  filtered_libraries_ch.single_sample
    .map{ it ->
      def processed = file("${params.results_dir}/${it.project_id}/${it.sample_id}/${it.unique_id}_processed.rds")
      def meta_json = file("${params.results_dir}/${it.project_id}/${it.sample_id}/${it.unique_id}_metadata.json")
      [it.unique_id, processed, meta_json ]
    }
    .subscribe{ unique_id, processed, meta_json ->
      if (!processed.exists() || !(processed.size() > 0)) {
        log.warn("Processed files do not exist for ${unique_id}. This library will not be included in the merged object.")
      }
      else if (!(meta_json.exists() && meta_json.size() > 0)) {
        log.warn("Metadata file does not exist for ${unique_id}. This library will not be included in the merged object.")
      }
      else if (getMetaVal(meta_json, "processed_cells") < 3) {
        log.warn("Library ${unique_id} has fewer than 3 cells. This library will not be included in the merged object.")
      }
    }

  present_libraries_ch = filtered_libraries_ch.single_sample
    .map{ it ->
      def processed = file("${params.results_dir}/${it.project_id}/${it.sample_id}/${it.unique_id}_processed.rds")
      def meta_json = file("${params.results_dir}/${it.project_id}/${it.sample_id}/${it.unique_id}_metadata.json")
      [it.project_id, it.unique_id, processed, meta_json]
    }
    // only include libraries that have been processed through scpca-nf and have at least 3 cells
    .filter{ _project_id, _unique_id, processed, meta_json ->
      (processed.exists() && processed.size() > 0 && getMetaVal(meta_json, "processed_cells") >= 3)
    }
    // remove metadata file
    .map{ project_id, unique_id, processed, _meta_json ->
      [project_id, unique_id, processed]
    }
    // only one row per library ID, this removes all the duplicates that may be present due to CITE/hashing
    .unique()
    // group tuple by project id: [project_id, [unique_id1, unique_id2, ...], [sce_file1, sce_file2, ...]]
    .groupTuple(by: 0)
    // remove any projects with more libraries present than the max allowed for merging
    .branch{ _project_id, unique_id_list, _sce_file_list ->
      oversized: unique_id_list.size() > params.max_merge_libraries
      within_limit: true
    }

  present_libraries_ch.oversized
    .subscribe{ project_id, unique_id_list, _sce_file_list ->
      log.warn("Not merging ${project_id} because it contains too many libraries (${unique_id_list.size()} > ${params.max_merge_libraries}).")
    }

  grouped_libraries_ch = present_libraries_ch.within_limit
    // add in boolean for if project contains samples with adt
    .map{ project_id, unique_id_list, sce_file_list ->
      [project_id, project_id in adt_projects.getVal(), unique_id_list, sce_file_list]
    }
    .branch{ it ->
      def merged_sce = file("${params.results_dir}/${it[0]}/merged/${it[0]}_merged.rds")
      has_merge: merged_sce.exists() && params.reuse_merge
      make_merge: true
    }

  pre_merged_ch = grouped_libraries_ch.has_merge
    .map{ project_id, has_adt, _unique_id_list, _sce_file_list ->
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

  // transpose so each h5ad file is checked separately
  anndata_format_ch = export_anndata.out
    .transpose()  // [merge_group_id, rna.h5ad], [merge_group_id, adt.h5ad]

  // check formatting of merged objects 
  check_sce(merged_ch, file(params.sce_format_reference_file))
  check_anndata(anndata_format_ch, file(params.anndata_format_reference_file))

  // collect all error files and concatenate to print to a formatting errors output file
  error_output_ch = check_sce.out
    .mix(check_anndata.out) // mix with the anndata error files

  error_files_ch = error_output_ch
    .collect{ merge_group_id, error_file -> error_file } // collect into a list of just the error files

  compile_errors(error_files_ch, params.merged_format_check_filename)

  // collect all error files and print out an error to the log file
  error_output_ch
    .filter{ meta, error_file -> error_file.size() > 0 } // only warn about objects with errors
    .subscribe{ meta, error_file -> 
      log.error "Formatting errors for ${error_file.text}"
    }
}
