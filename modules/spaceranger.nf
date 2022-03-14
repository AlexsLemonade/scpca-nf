#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process spaceranger{
  container params.SPACERANGER_CONTAINER
  publishDir "${params.outdir}/publish/${meta.project_id}/${meta.sample_id}"
  tag "${meta.run_id}-spatial" 
  label 'cpus_12'
  label 'mem_24'
  label 'disk_big'
  input:
    tuple val(meta), path(fastq_dir), file(image_file)
    path index
  output:
    tuple val(meta), path(spatial_out)
  script:
    spatial_out = "${meta.library_id}"
    out_id = "${meta.run_id}-spatial"
    meta.cellranger_index = index.fileName
    """
    spaceranger count \
      --id=${out_id} \
      --transcriptome=${index} \
      --fastqs=${fastq_dir} \
      --sample=${meta.cr_samples} \
      --localcores=${task.cpus} \
      --localmem=${task.memory.toGiga()} \
      --image=${image_file} \
      --slide=${meta.slide_serial_number} \
      --area=${meta.slide_section} 

    # make a new directory to hold only the outs file we want to publish 
    mkdir ${spatial_out}

    # move over needed files to outs directory 
    mv ${out_id}/outs/filtered_feature_bc_matrix ${spatial_out}
    mv ${out_id}/outs/raw_feature_bc_matrix ${spatial_out}
    mv ${out_id}/outs/spatial ${spatial_out}
    mv ${out_id}/outs/web_summary.html ${spatial_out}/${meta.library_id}_spaceranger_summary.html

    # move over versions file temporarily to be passed to metadata.json
    mv ${out_id}/_versions ${spatial_out}/spaceranger_versions.json
    mv ${out_id}/outs/metrics_summary.csv ${spatial_out}/spaceranger_metrics_summary.csv
    """
}

process spaceranger_metadata{
  container params.SCPCATOOLS_CONTAINER
  publishDir "${params.outdir}/publish/${meta.project_id}/${meta.sample_id}"
  input:
    tuple val(meta), path(spatial_out)
  output:
    tuple val(meta), path(metadata_json)
  script:
    metadata_json = "${meta.library_id}_metadata.json" 
    workflow_url = workflow.repository ?: workflow.manifest.homePage
    """
    generate_spaceranger_metadata.R \
      --library_id ${meta.library_id} \
      --sample_id ${meta.sample_id} \
      --unfiltered_barcodes_file "${spatial_out}/raw_feature_bc_matrix/barcodes.tsv.gz" \
      --filtered_barcodes_file "${spatial_out}/filtered_feature_bc_matrix/barcodes.tsv.gz" \
      --metrics_summary_file "${spatial_out}/spaceranger_metrics_summary.csv" \
      --spaceranger_versions_file "${spatial_out}/spaceranger_versions.json" \
      --metadata_json ${metadata_json} \
      --technology ${meta.technology} \
      --seq_unit ${meta.seq_unit} \
      --genome_assembly ${params.assembly} \
      --index_filename ${meta.cellranger_index} \
      --workflow_url "${workflow_url}" \
      --workflow_version "${workflow.revision}" \
      --workflow_commit "${workflow.commitId}"
    """
}

def getCRsamples(filelist){
  // takes a string with semicolon separated file names
  // returns just the 'sample info' portion of the file names,
  // as spaceranger would interpret them, comma separated
  fastq_files = filelist.tokenize(';').findAll{it.contains '.fastq.gz'}
  samples = []
  fastq_files.each{
    // append sample names to list, using regex to extract element before S001, etc.
    // [0] for the first match set, [1] for the first extracted element
    samples << (it =~ /^(.+)_S.+_L.+_[R|I].+.fastq.gz$/)[0][1]
  }
  // convert samples list to a `set` to remove duplicate entries,
  // then join to a comma separated string.
  return samples.toSet().join(',')
}


workflow spaceranger_quant{
    take: spatial_channel 
    // a channel with a map of metadata for each spatial library to process 
    main: 
        // create tuple of (metadata map, [])
        spaceranger_reads = spatial_channel
          // add sample names to metadata
          .map{it.cr_samples =  getCRsamples(it.files); it}
          // create tuple of [metadata, fastq dir, and path to image file]
          .map{meta -> tuple(meta,
                            file("${meta.files_directory}"),
                            file("${meta.files_directory}/*.jpg")
                            )}

        // run spaceranger
        spaceranger(spaceranger_reads, params.cellranger_index) \
          // generate metadata.json
          | spaceranger_metadata

    // tuple of metadata and path to spaceranger output directory
    emit: spaceranger.out
  
}
