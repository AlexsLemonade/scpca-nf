
// generate QC report from unfiltered and filtered SCE.rds files using scpcaTools

process sce_qc_report{
    container params.SCPCATOOLS_CONTAINER
    label 'mem_8'
    tag "${meta.library_id}"
    publishDir "${params.results_dir}/${meta.project_id}/${meta.sample_id}", mode: 'copy'
    input:
        tuple val(meta), path(unfiltered_rds), path(filtered_rds), path(processed_rds)
        tuple path(template_dir), val(template_file)
        tuple path(celltype_template_dir), val(celltype_template_file)
    output:
        tuple val(meta), path(unfiltered_out), path(filtered_out), path(processed_out), path(metadata_json), emit: data
        path qc_report, emit: report // TODO: CELLTYPE REPORT OUTPUT?
    script:
        qc_report = "${meta.library_id}_qc.html"
        template_path = "${template_dir}/${template_file}"
        metadata_json = "${meta.library_id}_metadata.json"
        workflow_url = workflow.repository ?: workflow.manifest.homePage
        workflow_version = workflow.revision ?: workflow.manifest.version
        
        // names for final output files
        unfiltered_out = "${meta.library_id}_unfiltered.rds"
        filtered_out = "${meta.library_id}_filtered.rds"
        processed_out = "${meta.library_id}_processed.rds"
        
        // check for cell types
        has_celltypes = params.perform_celltyping & (meta.submitter_cell_types_file | meta.singler_model_file | meta.cellassign_reference_file)
        celltype_report = "${meta.library_id}_celltype-report.html" // rendered HTML
        celltype_template_path = "${celltype_template_dir}/${celltype_template_file}" // template input

        """
        # move files for output
        if [ "${unfiltered_rds}" != "${unfiltered_out}" ]; then
            mv "${unfiltered_rds}" "${unfiltered_out}"
        fi
        if [ "${filtered_rds}" != "${filtered_out}" ]; then
            mv "${filtered_rds}" "${filtered_out}"
        fi
        if [ "${processed_rds}" != "${processed_out}" ]; then
            mv "${processed_rds}" "${processed_out}"
        fi

        # generate report and supplemental cell type report, if applicable
        sce_qc_report.R \
          --report_template "${template_path}" \
          --library_id "${meta.library_id}" \
          --sample_id "${meta.sample_id}" \
          --unfiltered_sce ${unfiltered_out} \
          --filtered_sce ${filtered_out} \
          --processed_sce ${processed_out} \
          --qc_report_file ${qc_report} \
          --celltypes_report_template "${celltype_template_path}" \
          ${has_celltypes ? "--cell_type_report_file ${celltype_report}" : ""} \
          --metadata_json ${metadata_json} \
          --technology "${meta.technology}" \
          --seq_unit "${meta.seq_unit}" \
          --genome_assembly "${meta.ref_assembly}" \
          --workflow_url "${workflow_url}" \
          --workflow_version "${workflow_version}" \
          --workflow_commit "${workflow.commitId}" \
          --seed "${params.seed}"
        """
    stub:
        unfiltered_out = "${meta.library_id}_unfiltered.rds"
        filtered_out = "${meta.library_id}_filtered.rds"
        processed_out = "${meta.library_id}_processed.rds"
        qc_report = "${meta.library_id}_qc.html"
        metadata_json = "${meta.library_id}_metadata.json"
        
        has_celltypes = params.perform_celltyping & (meta.submitter_cell_types_file | meta.singler_model_file | meta.cellassign_reference_file )
        celltype_report = "${meta.library_id}_celltype-report.html" // rendered HTML

        """
        touch ${unfiltered_out}
        touch ${filtered_out}
        touch ${processed_out}
        touch ${qc_report}
        echo '{"unfiltered_cells": 10, "filtered_cells": 10, "processed_cells": 10}' > ${metadata_json}
        
        if [ ${has_celltypes} == "true"]; 
          touch ${celltype_report}
        fi
        """
}
