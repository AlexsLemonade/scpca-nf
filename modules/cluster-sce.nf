// perform graph-based clustering on a processed SCE object
process cluster_sce{
    container params.SCPCATOOLS_CONTAINER
    label 'mem_8'
    tag "${meta.library_id}"
    input:
        tuple val(meta), path(unfiltered_rds), path(filtered_rds), path(processed_rds)
    output:
        tuple val(meta), path(unfiltered_rds), path(filtered_rds), path(clustered_rds)
    script:
        clustered_rds = "${meta.library_id}_clustered_processed.rds"
        """
        cluster_sce.R \
          --processed_sce_file ${processed_rds} \
          --cluster_algorithm ${params.cluster_algorithm} \
          --cluster_weighting ${params.cluster_weighting} \
          --nearest_neighbors ${params.nearest_neighbors} \
          ${params.seed ? "--random_seed ${params.seed}" : ""}
        mv ${processed_rds} ${clustered_rds}
        """
    stub:
        clustered_rds = "${meta.library_id}_clustered_processed.rds"
        """
        touch ${clustered_rds}
        """
}
