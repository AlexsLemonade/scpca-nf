// perform graph-based clustering on a processed SCE object
process cluster_sce {
  container params.SCPCATOOLS_SLIM_CONTAINER
  label 'mem_8'
  tag "${meta.library_id}"
  input:
    tuple val(meta), path(unfiltered_rds), path(filtered_rds), path(processed_rds)
  output:
    tuple val(meta), path(unfiltered_rds), path(filtered_rds), path(clustered_rds)
  script:
    clustered_rds = "${unique_id}_clustered.rds"
    """
    cluster_sce.R \
      --processed_sce_file ${processed_rds} \
      --output_sce_file ${clustered_rds} \
      --cluster_algorithm ${params.cluster_algorithm} \
      --cluster_weighting ${params.cluster_weighting} \
      --nearest_neighbors ${params.nearest_neighbors} \
      ${params.seed ? "--random_seed ${params.seed}" : ""}
    """
  stub:
    clustered_rds = "${unique_id}_clustered.rds"
    """
    touch ${clustered_rds}
    """
}
