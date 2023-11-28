
// Process to integrate SCEs with fastMNN
process integrate_fastmnn {
  container params.SCPCATOOLS_CONTAINER
  label 'mem_16'
  label 'cpus_4'
  input:
    tuple val(integration_group), path(merged_sce_file)
  output:
    tuple val(integration_group), path(integrated_sce_file)
  script:
    integrated_sce_file = "${integration_group}.rds"
    """
    integrate_sce.R \
      --input_sce_file "${merged_sce_file}" \
      --output_sce_file "${integrated_sce_file}" \
      --method "fastMNN" \
      --seed ${params.seed} \
      --threads ${task.cpus}
    """
  stub:
    integrated_sce_file = "${integration_group}.rds"
    """
    touch ${integrated_sce_file}
    """
}

// Process to integrate SCEs with Harmony
process integrate_harmony {
  container params.SCPCATOOLS_CONTAINER
  publishDir "${params.results_dir}/integration/${integration_group}"
  label 'mem_16'
  input:
    tuple val(integration_group), path(merged_sce_file)
  output:
    tuple val(integration_group), path(integrated_sce_file)
  script:
    integrated_sce_file = "${integration_group}.rds"
    """
    integrate_sce.R \
      --input_sce_file "${merged_sce_file}" \
      --output_sce_file "${integrated_sce_file}" \
      --method "harmony" \
      --seed ${params.seed}
    """
  stub:
    integrated_sce_file = "${integration_group}.rds"
    """
    touch ${integrated_sce_file}
    """
}
