#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process spaceranger{
  container params.SPACERANGER_CONTAINER
  publishDir "${meta.spaceranger_publish_dir}"
  tag "${meta.run_id}-spatial" 
  label 'cpus_12'
  label 'mem_24'
  label 'disk_big'
  input:
    tuple val(meta), path(fastq_dir), file(image_file)
    path index
  output:
    tuple val(meta), path(out_id)
  script:
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

    # remove bam and bai files
    rm ${out_id}/outs/*.bam*
    """
}

process spaceranger_publish{
  container params.SCPCATOOLS_CONTAINER
  publishDir "${params.outdir}/publish/${meta.project_id}/${meta.sample_id}"
  input:
    tuple val(meta), path(spatial_out)
  output:
    tuple val(meta), path(spatial_publish_dir), path(metadata_json)
  script:
    spatial_publish_dir = "${meta.library_id}_spatial"
    metadata_json = "${meta.library_id}_metadata.json" 
    workflow_url = workflow.repository ?: workflow.manifest.homePage
    """
    # make a new directory to hold only the outs file we want to publish 
    mkdir ${spatial_publish_dir}

    # move over needed files to outs directory 
    mv ${spatial_out}/outs/filtered_feature_bc_matrix ${spatial_publish_dir}
    mv ${spatial_out}/outs/raw_feature_bc_matrix ${spatial_publish_dir}
    mv ${spatial_out}/outs/spatial ${spatial_publish_dir}
    mv ${spatial_out}/outs/web_summary.html ${spatial_publish_dir}/${meta.library_id}_spaceranger_summary.html

    generate_spaceranger_metadata.R \
      --library_id ${meta.library_id} \
      --sample_id ${meta.sample_id} \
      --unfiltered_barcodes_file "${spatial_publish_dir}/raw_feature_bc_matrix/barcodes.tsv.gz" \
      --filtered_barcodes_file "${spatial_publish_dir}/filtered_feature_bc_matrix/barcodes.tsv.gz" \
      --metrics_summary_file "${spatial_out}/outs/metrics_summary.csv" \
      --spaceranger_versions_file "${spatial_out}/_versions" \
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

def getCRsamples(files_dir){
  // takes a string with semicolon separated file names
  // returns just the 'sample info' portion of the file names,
  // as spaceranger would interpret them, comma separated
  fastq_files = file(files_dir).list().findAll{it.contains '.fastq.gz'}
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
        //grouped_spatial_channel = spatial_channel
          // add files, sample names and spatial output directory to metadata
          //.map{meta -> tuple(meta,
          /*                    meta.library_id,
                             file("${meta.files_directory}/*_R*_*.fastq.gz")
                             )}
          // one fastq file per line 
          .transpose()
          // add new entry with 
          .map{it.cr_samples = (it =~ /^(.+)_S.+_L.+_[R|I].+.fastq.gz$/)[0][1];
               it}
          .groupTuple(by:1) // group by library ID 
          .map{[it[0], it[1], it[2]]} // remove files and put into one tuple  */

          
        spatial_channel = spatial_channel
          .map{it.cr_samples = getCRsamples(it.files_directory);
               it.spaceranger_publish_dir =  "${params.outdir}/internal/spaceranger/${it.library_id}";
               it.spaceranger_results_dir = "${it.spaceranger_publish_dir}/${it.run_id}-spatial";
               it}
          .branch{
            has_spatial: !params.repeat_mapping & file(it.spaceranger_results_dir).exists()
            make_spatial: true
           }

          // create tuple of [metadata, fastq dir, and path to image file]
        spaceranger_reads = spatial_channel.make_spatial
          .map{meta -> tuple(meta,
                            file("${meta.files_directory}"),
                            file("${meta.files_directory}/*.jpg")
                            )}

        // run spaceranger
        spaceranger(spaceranger_reads, params.cellranger_index)

        // gather spaceranger output for completed libraries
        spaceranger_quants_ch = spatial_channel.has_spatial
          .map{meta -> tuple(meta,
                             file("${meta.spaceranger_results_dir}")
                             )}

        grouped_spaceranger_ch = spaceranger.out.mix(spaceranger_quants_ch)

          // generate metadata.json
        spaceranger_publish(grouped_spaceranger_ch)

    // tuple of metadata, path to spaceranger output directory, and path to metadata json file
    emit: spaceranger_publish.out
  
}
