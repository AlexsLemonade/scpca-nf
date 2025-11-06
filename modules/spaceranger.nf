#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process spaceranger {
  container params.SPACERANGER_CONTAINER
  publishDir "${meta.spaceranger_publish_dir}", mode: 'copy'
  tag "${meta.run_id}-spatial"
  label 'cpus_12'
  label 'mem_24'
  label 'disk_big'
  input:
    tuple val(meta), path(fastq_dir), path(image_file), path(index)
  output:
    tuple val(meta), path(out_id)
  script:
    out_id = file(meta.spaceranger_results_dir).name
    meta += Utils.getVersions(workflow, nextflow)
    meta_json = Utils.makeJson(meta)
    """
    spaceranger count \
      --id=${out_id} \
      --transcriptome=${index} \
      --fastqs=${fastq_dir} \
      --sample=${meta.cr_samples} \
      --localcores=${task.cpus} \
      --localmem=${task.memory.toGiga()} \
      --image=${image_file} \
      --slide=${meta.slide_serial_number ?: "NA"} \
      --area=${meta.slide_section ?: "NA"}

    # write metadata
    echo '${meta_json}' > ${out_id}/scpca-meta.json

    # remove Space Ranger intermediates directory
    rm -rf ${out_id}/SPATIAL_RNA_COUNTER_CS

    # remove bam and bai files
    rm ${out_id}/outs/*.bam*
    """
  stub:
    out_id = file(meta.spaceranger_results_dir).name
    meta += Utils.getVersions(workflow, nextflow)
    meta_json = Utils.makeJson(meta)
    """
    mkdir -p ${out_id}/outs
    echo '${meta_json}' > ${out_id}/scpca-meta.json
    """
}

process spaceranger_publish {
  container params.SCPCATOOLS_SLIM_CONTAINER
  tag "${meta.library_id}"
  publishDir "${params.results_dir}/${meta.project_id}/${meta.sample_id}", mode: 'copy'
  input:
    tuple val(meta), path(spatial_out)
  output:
    tuple val(meta), path(spatial_publish_dir), path(metadata_json)
  script:
    spatial_publish_dir = "${meta.library_id}_spatial"
    metadata_json = "${spatial_publish_dir}/${meta.library_id}_metadata.json"
    workflow_url = workflow.repository ?: workflow.manifest.homePage
    cellranger_index_name = file(meta.cellranger_index).name
    """
    # make a new directory to hold only the outs file we want to publish
    mkdir ${spatial_publish_dir}

    # move over needed files to outs directory
    mv ${spatial_out}/outs/filtered_feature_bc_matrix ${spatial_publish_dir}
    mv ${spatial_out}/outs/raw_feature_bc_matrix ${spatial_publish_dir}
    mv ${spatial_out}/outs/spatial ${spatial_publish_dir}
    mv ${spatial_out}/outs/web_summary.html ${spatial_publish_dir}/${meta.library_id}_spaceranger-summary.html

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
      --genome_assembly ${meta.ref_assembly} \
      --index_filename ${cellranger_index_name} \
      --workflow_url "${workflow_url}" \
      --workflow_version "${workflow.revision}" \
      --workflow_commit "${workflow.commitId}"
    """
  stub:
    spatial_publish_dir = "${meta.library_id}_spatial"
    metadata_json = "${spatial_publish_dir}/${meta.library_id}_metadata.json"
    """
    mkdir -p ${spatial_publish_dir}
    echo '{}' > ${metadata_json}
    """
}

def getCRsamples(files_dir) {
  // takes the path to the directory holding the fastq files for each sample
  // returns just the 'sample info' portion of the file names,
  // as spaceranger would interpret them, comma separated
  if (!files_dir) { // return empty string if no files directory
    return ""
  }
  def fastq_files = files("${files_dir}/*.fastq.gz")
  def samples = []
  fastq_files.each{ it ->
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
    spatial_channel = spatial_channel
      // add sample names and spatial output directory to metadata
      .map{ meta_in ->
        def meta = meta_in.clone()
        meta.cr_samples = getCRsamples(meta.files_directory)
        meta.spaceranger_publish_dir =  "${params.checkpoints_dir}/spaceranger/${meta.library_id}"
        meta.spaceranger_results_dir = "${meta.spaceranger_publish_dir}/${meta.run_id}-spatial"
        meta // return modified meta object
      }
      .branch{ it ->
        def stored_ref_assembly = Utils.getMetaVal(file("${it.spaceranger_results_dir}/scpca-meta.json"), "ref_assembly")
        def stored_tech = Utils.getMetaVal(file("${it.spaceranger_results_dir}/scpca-meta.json"), "technology").toLowerCase()
        make_spatial: (
          // input files exist
          it.files_directory && file(it.files_directory, type: "dir").exists() && (
            params.repeat_mapping
            || !file(it.spaceranger_results_dir).exists()
            // or the technology has changed (to ensure re-mapping if tech was updated)
            || it.technology.toLowerCase() != stored_tech
            || it.ref_assembly != stored_ref_assembly
          )
        )
        has_spatial: file(it.spaceranger_results_dir).exists()
        missing_inputs: true
      }

    // send run ids in spatial_channel.missing_inputs to log
    spatial_channel.missing_inputs
      .subscribe{ it ->
        log.error("The expected input files or Space Ranger results files for ${it.run_id} are missing.")
      }

    // create tuple of [metadata, fastq dir, and path to image file]
    spaceranger_reads = spatial_channel.make_spatial
      .map{ meta ->
        [
          meta,
          file(meta.files_directory, type: 'dir'),
          file("${meta.files_directory}/*.jpg"),
          file(meta.cellranger_index, type: 'dir')
        ]
    }

    // run spaceranger
    spaceranger(spaceranger_reads)

    // gather spaceranger output for completed libraries
    // make a list of metadata (read from prior output) and prior results directory
    spaceranger_quants_ch = spatial_channel.has_spatial
      .map{ meta ->
        [
          Utils.readMeta(file("${meta.spaceranger_results_dir}/scpca-meta.json")),
          file(meta.spaceranger_results_dir, type: 'dir')
        ]
      }

    grouped_spaceranger_ch = spaceranger.out.mix(spaceranger_quants_ch)

    // generate metadata.json
    spaceranger_publish(grouped_spaceranger_ch)

  // tuple of metadata, path to spaceranger output directory, and path to metadata json file
  emit: spaceranger_publish.out

}
