
// generate QC report from SCE files and publish SCE files and JSONs

include { pullthroughContainer } from '../lib/utils.nf'

process qc_publish_sce {
  container "${pullthroughContainer(params.scpcatools_reports_container, params.pullthrough_registry)}"
  label 'mem_16'
  tag "${meta.unique_id}"
  publishDir "${params.results_dir}/${meta.project_id}/${meta.sample_id}", mode: 'copy'
  input:
    tuple val(meta), path(unfiltered_rds), path(filtered_rds), path(processed_rds), path(infercnv_heatmap_file)
    tuple path(template_dir), val(qc_template_file), val(celltype_template_file)
    val(perform_celltyping)
    path(validation_groups_file)
    path(validation_markers_file)
    path(validation_palette_file)
  output:
    tuple val(meta), path(unfiltered_out), path(filtered_out), path(processed_out), path(metadata_json), emit: data
    path qc_report, emit: report
    path celltype_report, emit: celltype_report, optional: true
    path metrics_json, emit: metrics
  script:
    workflow_url = workflow.repository ?: workflow.manifest.homePage
    workflow_version = workflow.revision ?: workflow.manifest.version

    // determine if we have a usable heatmap file
    has_infercnv = infercnv_heatmap_file && infercnv_heatmap_file.isFile() && infercnv_heatmap_file.size() > 0

    // names for final output files
    unfiltered_out = "${meta.unique_id}_unfiltered.rds"
    filtered_out = "${meta.unique_id}_filtered.rds"
    processed_out = "${meta.unique_id}_processed.rds"
    qc_report = "${meta.unique_id}_qc.html"
    metadata_json = "${meta.unique_id}_metadata.json"
    metrics_json = "${meta.unique_id}_metrics.json"

    // check for cell types
    // only provide report template if cell typing was performed and either singler or cellassign was used
    // note that we use the value perform_celltyping, not the param here
    has_celltypes = perform_celltyping && (meta.singler_model_file || meta.cellassign_reference_file)
    celltype_report = "${meta.unique_id}_celltype-report.html" // rendered HTML

    """
    # copy files for output if needed
    if [ "${unfiltered_rds}" != "${unfiltered_out}" ]; then
      cp "${unfiltered_rds}" "${unfiltered_out}"
    fi
    if [ "${filtered_rds}" != "${filtered_out}" ]; then
      cp "${filtered_rds}" "${filtered_out}"
    fi
    if [ "${processed_rds}" != "${processed_out}" ]; then
      cp "${processed_rds}" "${processed_out}"
    fi

    # generate report and supplemental cell type report, if applicable
    sce_qc_report.R \
      --report_template "${template_dir / qc_template_file}" \
      --validation_groups_file ${validation_groups_file} \
      --validation_markers_file ${validation_markers_file} \
      --validation_palette_file ${validation_palette_file} \
      ${has_infercnv ? "--infercnv_heatmap_file ${infercnv_heatmap_file}" : ""} \
      --library_id "${meta.library_id}" \
      --sample_id "${meta.sample_id}" \
      --project_id "${meta.project_id}" \
      --unfiltered_sce ${unfiltered_out} \
      --filtered_sce ${filtered_out} \
      --processed_sce ${processed_out} \
      --qc_report_file ${qc_report} \
      --celltype_report_template "${template_dir / celltype_template_file}" \
      ${has_celltypes ? "--celltype_report_file ${celltype_report}" : ""} \
      --metadata_json ${metadata_json} \
      --technology "${meta.technology}" \
      --seq_unit "${meta.seq_unit}" \
      --genome_assembly "${meta.ref_assembly}" \
      --infercnv_min_reference_cells ${params.infercnv_min_reference_cells} \
      --workflow_url "${workflow_url}" \
      --workflow_version "${workflow_version}" \
      --workflow_commit "${workflow.commitId}" \
      --seed "${params.seed}"

    sce_metrics.R \
      --metadata_json "${metadata_json}" \
      --unfiltered_sce "${unfiltered_out}" \
      --filtered_sce "${filtered_out}" \
      --processed_sce "${processed_out}" \
      --metrics_json "${metrics_json}"
    """
  stub:
    unfiltered_out = "${meta.unique_id}_unfiltered.rds"
    filtered_out = "${meta.unique_id}_filtered.rds"
    processed_out = "${meta.unique_id}_processed.rds"
    qc_report = "${meta.unique_id}_qc.html"
    metadata_json = "${meta.unique_id}_metadata.json"
    metrics_json = "${meta.unique_id}_metrics.json"

    has_celltypes = params.perform_celltyping && (meta.singler_model_file || meta.cellassign_reference_file)
    celltype_report = "${meta.unique_id}_celltype-report.html" // rendered HTML

    """
    touch ${unfiltered_out}
    touch ${filtered_out}
    touch ${processed_out}
    touch ${qc_report}
    touch ${metrics_json}
    ${has_celltypes ? "touch ${celltype_report}" : ""}

    echo '{"unfiltered_cells": 10, "filtered_cells": 10, "processed_cells": 10}' > ${metadata_json}
    """
}
