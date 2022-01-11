// generate unfiltered and filtered RDS files using scpcaTools

// RNA only libraries
process make_unfiltered_sce{
    container params.SCPCATOOLS_CONTAINER
    label 'mem_8'
    publishDir "${params.outdir}/publish/${meta.project_id}/${meta.sample_id}"
    input: 
        tuple val(meta), path(alevin_dir)
        path(mito)
        path(gtf)
    output:
        tuple val(meta), path(unfiltered_rds)
    script:
        unfiltered_rds = "${meta.library_id}_unfiltered.rds"
        """
        generate_unfiltered_sce.R \
          --seq_unit ${meta.seq_unit} \
          --alevin_dir ${alevin_dir} \
          --unfiltered_file ${unfiltered_rds} \
          --mito_file ${mito} \
          --gtf_file ${gtf} \
          --technology ${meta.technology}
        """
}

// channels with RNA and feature data
process make_merged_unfiltered_sce{
    label 'mem_8'
    container params.SCPCATOOLS_CONTAINER
    publishDir "${params.outdir}/publish/${meta.project_id}/${meta.sample_id}"
    input: 
        tuple val(feature_meta), path(feature_alevin_dir), val (meta), path(alevin_dir)
        path(mito)
        path(gtf)
    output:
        tuple val(meta), path(unfiltered_rds)
    script:
        unfiltered_rds = "${meta.library_id}_unfiltered.rds"
        // add feature metadata as an element of the main meta object
        meta['feature_type'] = feature_meta.technology.split('_')[0]
        meta['feature_meta'] = feature_meta
        
        """
        generate_unfiltered_sce.R \
          --seq_unit ${meta.seq_unit} \
          --alevin_dir ${alevin_dir} \
          --feature_dir ${feature_alevin_dir} \
          --feature_name ${meta.feature_type} \
          --unfiltered_file ${unfiltered_rds} \
          --mito_file ${mito} \
          --gtf_file ${gtf} \
          --technology ${meta.technology}
        """
}

process filter_sce{
    container params.SCPCATOOLS_CONTAINER
    label 'mem_8'
    publishDir "${params.outdir}/publish/${meta.project_id}/${meta.sample_id}"
    input: 
        tuple val(meta), path(unfiltered_rds)
    output:
        tuple val(meta), path(unfiltered_rds), path(filtered_rds)
    script:
        filtered_rds = "${meta.library_id}_filtered.rds"
        """
        filter_sce_rds.R \
          --unfiltered_file ${unfiltered_rds} \
          --filtered_file ${filtered_rds} \
          --lower ${params.emptydrops_lower} \
          ${params.seed ? "--random_seed ${params.seed}" : ""}
        """
}

workflow generate_sce {
  // generate rds files for RNA-only samples
  take: quant_channel
  main:
    make_unfiltered_sce(quant_channel, params.mito_file, params.ref_gtf) \
      | filter_sce

  emit: filter_sce.out
  // a tuple of meta and the filtered and unfiltered rds files
}

workflow generate_merged_sce {
  // generate rds files for feature + quant samples
  // input is a channel with feature_meta, feature_quantdir, rna_meta, rna_quantdir
  take: feature_quant_channel
  main:
    make_merged_unfiltered_sce(feature_quant_channel, params.mito_file, params.ref_gtf) \
      | filter_sce

  emit: filter_sce.out
  // a tuple of meta and the filtered and unfiltered rds files
}
