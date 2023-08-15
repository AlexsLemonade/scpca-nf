// perform graph-based clustering on a processed SCE object
// this process also does the RDS file export to the publishDir


process cluster_sce{
    container params.SCPCATOOLS_CONTAINER
    label 'mem_8'
    tag "${meta.library_id}"
    publishDir "${params.results_dir}/${meta.project_id}/${meta.sample_id}", mode: 'copy'
    input:
        tuple val(meta), path(unfiltered_rds), path(filtered_rds), path(processed_rds)
    output:
        tuple val(meta), path(unfiltered_rds), path(filtered_rds), path(processed_rds)
    script:
        """
        cluster_sce.R \
          --processed_sce_file ${processed_rds} \
          --cluster_algorithm ${params.cluster_algorithm} \
          --cluster_weighting ${params.cluster_weighting} \
          --nearest_neighbors ${params.nearest_neighbors} \
          ${params.seed ? "--random_seed ${params.seed}" : ""}
        """
    stub:
        processed_rds = "${meta.library_id}_processed.rds"
        """
        touch ${processed_rds}
        """
}
