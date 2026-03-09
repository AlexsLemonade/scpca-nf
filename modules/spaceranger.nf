
include { parseNA; getVersions; makeJson; readMeta; getMetaVal; pullthroughContainer; getImageFiles } from '../lib/utils.nf'

process spaceranger {
  container "${pullthroughContainer(params.spaceranger_container, params.pullthrough_registry)}"
  publishDir "${meta.spaceranger_publish_dir}", mode: 'copy'
  tag "${meta.run_id}-spatial"
  label 'cpus_12'
  label 'mem_32'
  label 'disk_big'
  input: 
    tuple val(meta), path(index), path(probeset_file),
      path(fastq_dir), path(cytaimage_file), path(image_file), val(image_type)
  output:
    tuple val(meta), path(out_id)
  script:
    out_id = file(meta.spaceranger_results_dir).name
    meta += getVersions()
    meta_json = makeJson(meta)

    // TODO: We are currently assuming that only one of the non-cytassist images should exist, 
    //  and therefore we only want to provide a single argument
    image_arg = image_type == "image"           ? "--image ${image_file}" :
                image_type == "colorizedimage"  ? "--colorizedimage ${image_file}" :
                image_type == "darkimage"       ? "--darkimage ${[image_file].flatten().join(',')}" :
                ""

    """
    spaceranger count \
      --id=${out_id} \
      --transcriptome=${index} \
      --fastqs=${fastq_dir} \
      --sample=${meta.cr_samples} \
      --localcores=${task.cpus} \
      --localmem=${task.memory.toGiga()} \
      --slide=${meta.slide_serial_number} \
      --area=${meta.slide_section} \
      --create-bam false \
      ${probeset_file ? "--probe-set ${probeset_file}" : ""} \
      ${cytaimage_file ? "--cytaimage ${cytaimage_file}" : ""} \
      ${image_arg}

    # write metadata
    echo '${meta_json}' > ${out_id}/scpca-meta.json

    # TODO: Why is this still present in the checkpoints directory???
    # TODO: Is there more we'd like to remove?
    # remove Space Ranger intermediates directory
    rm -rf ${out_id}/SPATIAL_RNA_COUNTER_CS
    """
  stub:
    out_id = file(meta.spaceranger_results_dir).name
    meta += getVersions()
    meta_json = makeJson(meta)
    """
    mkdir -p ${out_id}/outs
    echo '${meta_json}' > ${out_id}/scpca-meta.json
    """
}

// TODO: need to accommodate HD outputs which are slightly different
process spaceranger_publish {
  container "${pullthroughContainer(params.scpcatools_slim_container, params.pullthrough_registry)}"
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

    # copy over needed files to outs directory
    cp -r ${spatial_out}/outs/filtered_feature_bc_matrix ${spatial_publish_dir}
    cp -r ${spatial_out}/outs/raw_feature_bc_matrix ${spatial_publish_dir}
    cp -r ${spatial_out}/outs/spatial ${spatial_publish_dir}
    cp ${spatial_out}/outs/web_summary.html ${spatial_publish_dir}/${meta.library_id}_spaceranger-summary.html

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
    // use .name to ensure we aren't getting a full path
    samples << (it.name =~ /^(.+)_S.+_L.+_[R|I].+.fastq.gz$/)[0][1]
  }
  // convert samples list to a `set` to remove duplicate entries,
  // then join to a comma separated string.
  return samples.toSet().join(',')
}


workflow spaceranger_quant{
  take: spatial_channel // a channel with a map of metadata for each spatial library to process
  main:
    // techs that use the probe file
    def cytassist_probe_techs = ['visium2', 'visium_hd_3prime']

    spatial_channel = spatial_channel
      // add sample names and spatial output directory to metadata
      .map{ meta_in ->
        def meta = meta_in.clone()
        meta.cr_samples = getCRsamples("${meta.files_directory}/fastq")
        meta.spaceranger_publish_dir =  "${params.checkpoints_dir}/spaceranger/${meta.library_id}"
        meta.spaceranger_results_dir = "${meta.spaceranger_publish_dir}/${meta.run_id}-spatial"

        meta // return modified meta object
      }
      .branch{ it ->
        def stored_ref_assembly = getMetaVal(file("${it.spaceranger_results_dir}/scpca-meta.json"), "ref_assembly")
        def stored_tech = getMetaVal(file("${it.spaceranger_results_dir}/scpca-meta.json"), "technology") ?: ""
        // branch for invalid cases
        missing_slide_serial: !parseNA(it.slide_serial_number)
        missing_slide_section: !parseNA(it.slide_section)
        make_spatial: (
          // input files exist
          it.files_directory && file(it.files_directory, type: "dir").exists() && (
            params.repeat_mapping
            || !file(it.spaceranger_results_dir).exists()
            // or the technology has changed (to ensure re-mapping if tech was updated)
            || it.technology.toLowerCase() != stored_tech.toLowerCase()
            || it.ref_assembly != stored_ref_assembly
          )
        )
        has_spatial: file(it.spaceranger_results_dir).exists()
        missing_inputs: true
      }

    // TODO: Should these be warn or error?
    // warn about missing spatial information
    spatial_channel.missing_slide_serial.subscribe{ it ->
      log.warn("Run '${it.run_id}' is missing a slide serial number and will not be processed.")
    }
    spatial_channel.missing_slide_section.subscribe{ it ->
      log.warn("Run '${it.run_id}' is missing a slide section (area) and will not be processed.")
    }
    spatial_channel.missing_inputs.subscribe{ it ->
      log.error("The expected input files or Space Ranger results files for ${it.run_id} are missing.")
    }
  

    // create tuple of [metadata, index, probeset file, fastq_dir, cytaimage file, image file, image type]
    spaceranger_reads = spatial_channel.make_spatial
      .map{ meta ->
        // probeset logic
        def use_probeset = meta.technology in cytassist_probe_techs
        def probeset_file = use_probeset ? file(meta.cytassist_probe) : []

        // image logic: first three must be a single file, but darkimage can be several
        def cytaimage_file = getImageFiles("${meta.files_directory}/cytaimage")
        def image_file = getImageFiles("${meta.files_directory}/image")
        def colorizedimage_file = getImageFiles("${meta.files_directory}/colorizedimage")
        def darkimage_files = getImageFiles("${meta.files_directory}/darkimage/*", 6) // max 6 files: https://www.10xgenomics.com/support/software/space-ranger/latest/analysis/inputs/image-image-recommendation

        // TODO: We are currently assuming that only one of the non-cytassist images should exist
        def image_type = image_file ? "image" :
                         colorizedimage_file ? "colorizedimage" :
                         darkimage_files ? "darkimage" :
                         ""

        [
          meta,
          file(meta.cellranger_index, type: 'dir'),
          probeset_file,
          file("${meta.files_directory}/fastq", type: 'dir'),
          cytaimage_file, 
          image_file, 
          image_type   
        ]
    }

    // run spaceranger
    spaceranger(spaceranger_reads)

    // gather spaceranger output for completed libraries
    // make a list of metadata (read from prior output) and prior results directory
    spaceranger_quants_ch = spatial_channel.has_spatial
      .map{ meta ->
        [
          readMeta(file("${meta.spaceranger_results_dir}/scpca-meta.json")),
          file(meta.spaceranger_results_dir, type: 'dir')
        ]
      }

    grouped_spaceranger_ch = spaceranger.out.mix(spaceranger_quants_ch)

    // generate metadata.json
    spaceranger_publish(grouped_spaceranger_ch)

  // tuple of metadata, path to spaceranger output directory, and path to metadata json file
  emit: spaceranger_publish.out

}
