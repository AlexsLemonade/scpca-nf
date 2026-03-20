
include { getVersions; makeJson; readMeta; getMetaVal; pullthroughContainer; getImageFiles } from '../lib/utils.nf'

process spaceranger {
  container "${pullthroughContainer(params.spaceranger_container, params.pullthrough_registry)}"
  publishDir "${meta.spaceranger_checkpoint_dir}", mode: 'copy'
  tag "${meta.run_id}-spatial"
  label 'cpus_12'
  label 'mem_32'
  label 'disk_big'
  input: 
    tuple val(meta), path(index), path(probeset_file),
      path(fastq_dir), path(cytaimage_file), path(image_file)
  output:
    tuple val(meta), path(out_id)
  script:
    out_id = file(meta.spaceranger_results_dir).name
    meta += getVersions()
    meta_json = makeJson(meta)
    image_arg = meta.visium_image_type ? "--${meta.visium_image_type} ${image_file.join(',')}" : "" // join in case of multiple darkimages
    
    // may help avoid OOM errors, but needs to be a rounded integer sans decimal
    spaceranger_mem = Math.ceil(task.memory.toGiga() * 0.9) as int
    """
    spaceranger count \
      --id=${out_id} \
      --transcriptome=${index} \
      --fastqs=${fastq_dir} \
      --sample=${meta.cr_samples} \
      --localcores=${task.cpus} \
      --localmem=${spaceranger_mem} \
      --slide=${meta.slide_serial_number} \
      --area=${meta.slide_section} \
      --create-bam=false \
      ${probeset_file ? "--probe-set ${probeset_file}" : ""} \
      ${cytaimage_file ? "--cytaimage ${cytaimage_file}" : ""} \
      ${image_arg}

    # remove Space Ranger intermediates directory 
    # do this before saving the json
    rm -rf ${out_id}/SPATIAL_RNA_COUNTER_CS
    
    # write metadata
    echo '${meta_json}' > ${out_id}/scpca-meta.json
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


// this process is identical to spaceranger, but with more memory for HD & HD 3' datasets
// use 12 CPU/64 RAM if no brightfield image, and 16 CPU/128 RAM if there is to accommodate segmentation memory needs
process spaceranger_hd {
  container "${pullthroughContainer(params.spaceranger_container, params.pullthrough_registry)}"
  publishDir "${meta.spaceranger_checkpoint_dir}", mode: 'copy'
  tag "${meta.run_id}-spatial"
  memory {
    def mem = meta.visium_image_type == "image" ? 64.GB : 32.GB 
    mem + mem * task.attempt
  }
  cpus {
    meta.visium_image_type == "image" ? 16 : 12
  }
  label 'disk_big'
  input: 
    tuple val(meta), path(index), path(probeset_file),
      path(fastq_dir), path(cytaimage_file), path(image_file)
  output:
    tuple val(meta), path(out_id)
  script:
    out_id = file(meta.spaceranger_results_dir).name
    meta += getVersions()
    meta_json = makeJson(meta)
    image_arg = meta.visium_image_type ? "--${meta.visium_image_type} ${image_file.join(',')}" : "" // join in case of multiple darkimages
    
    // may help avoid OOM errors, but needs to be a rounded integer sans decimal
    spaceranger_mem = Math.ceil(task.memory.toGiga() * 0.9) as int
    """
    spaceranger count \
      --id=${out_id}-uncompressed \
      --transcriptome=${index} \
      --fastqs=${fastq_dir} \
      --sample=${meta.cr_samples} \
      --localcores=${task.cpus} \
      --localmem=${spaceranger_mem} \
      --slide=${meta.slide_serial_number} \
      --area=${meta.slide_section} \
      --create-bam=false \
      ${probeset_file ? "--probe-set ${probeset_file}" : ""} \
      ${cytaimage_file ? "--cytaimage ${cytaimage_file}" : ""} \
      ${image_arg}

    mkdir -p ${out_id}/outs 
  
    # compress outs/ into out_id
    tar -czf ${out_id}/${meta.library_id}_spatial.tar.gz ${out_id}-uncompressed/outs

    # add indicator if the segmented directory exists
    if [ -d "${out_id}-uncompressed/outs/segmented_outputs" ]; then
      touch ${out_id}/segmented_exists.txt
    fi

    # get rest of the things into out_id
    cp ${out_id}-uncompressed/outs/binned_outputs/square_008um/raw_feature_bc_matrix/barcodes.tsv.gz ${out_id}/raw-barcodes.tsv.gz
    cp ${out_id}-uncompressed/outs/binned_outputs/square_008um/filtered_feature_bc_matrix/barcodes.tsv.gz ${out_id}/filtered-barcodes.tsv.gz
    cp ${out_id}-uncompressed/outs/metrics_summary.csv ${out_id}/outs/metrics_summary.csv
    cp ${out_id}-uncompressed/outs/web_summary.html ${out_id}/outs/web_summary.html
    cp ${out_id}-uncompressed/_versions ${out_id}/_versions

    # write metadata
    echo '${meta_json}' > ${out_id}/scpca-meta.json
    """
  stub:
    out_id = file(meta.spaceranger_results_dir).name
    meta += getVersions()
    meta_json = makeJson(meta)
    """
    mkdir -p ${out_id}/outs
    echo '${meta_json}' > ${out_id}/scpca-meta.json
    if [ "${meta.visium_image_type}" == "image" ]; then
      touch ${out_id}/segmented_exists.txt
    fi
    """
}

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

    is_hd = meta.technology.contains("visium-hd")
    """
    # make a new directory to hold only the outs file we want to publish
    mkdir ${spatial_publish_dir}

    # copy over summary html, present for all technologies
    cp ${spatial_out}/outs/web_summary.html ${spatial_publish_dir}/${meta.library_id}_spaceranger-summary.html

    # copy over files that depend on the technology, and define associated inputs to the R script
    if [[ ${is_hd} == "true" ]]; then
       cp ${spatial_out}/${meta.library_id}_spatial.tar.gz ${spatial_publish_dir}/

      unfiltered_barcodes_file="${spatial_out}/raw-barcodes.tsv.gz"
      filtered_barcodes_file="${spatial_out}/filtered-barcodes.tsv.gz"
    else 
      cp -r ${spatial_out}/outs/raw_feature_bc_matrix ${spatial_publish_dir}
      cp -r ${spatial_out}/outs/filtered_feature_bc_matrix ${spatial_publish_dir}
      cp -r ${spatial_out}/outs/spatial ${spatial_publish_dir}
      
      unfiltered_barcodes_file="${spatial_out}/outs/raw_feature_bc_matrix/barcodes.tsv.gz"
      filtered_barcodes_file="${spatial_out}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"
    fi

    generate_spaceranger_metadata.R \
      --library_id ${meta.library_id} \
      --sample_id ${meta.sample_id} \
      --unfiltered_barcodes_file \${unfiltered_barcodes_file} \
      --filtered_barcodes_file \${filtered_barcodes_file} \
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
  take: 
    spatial_channel // a channel with a map of metadata for each spatial library to process
    spatial_techs   // list of possible spatial technologies
  main:
    // techs that require a probe file
    def visium_probe_techs = spatial_techs.findAll{ it =~ /_v\d\.*\d*$/}

    // for image handling, techs that are _not_ cytassist
    def non_cytassist_techs = ['visium', 'visium1', 'visium1_v1']

    spatial_channel = spatial_channel
      // add sample names and spatial output directory to metadata
      .map{ meta_in ->
        def meta = meta_in.clone()
        meta.cr_samples = getCRsamples("${meta.files_directory}/fastq")
        meta.spaceranger_checkpoint_dir =  "${params.checkpoints_dir}/spaceranger/${meta.library_id}"
        meta.spaceranger_results_dir = "${meta.spaceranger_checkpoint_dir}/${meta.run_id}-spatial"

        meta // return modified meta object
      }
      .branch{ it ->
        def stored_ref_assembly = getMetaVal(file("${it.spaceranger_results_dir}/scpca-meta.json"), "ref_assembly")
        def stored_tech = getMetaVal(file("${it.spaceranger_results_dir}/scpca-meta.json"), "technology") ?: ""
        // branch for invalid cases
        missing_slide_serial: !it.slide_serial_number
        missing_slide_section: !it.slide_section
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

    spatial_channel.missing_slide_serial.subscribe{ it ->
      log.error("Run '${it.run_id}' is missing a slide serial number and will not be processed.")
    }
    spatial_channel.missing_slide_section.subscribe{ it ->
      log.error("Run '${it.run_id}' is missing a slide section (area) and will not be processed.")
    }
    spatial_channel.missing_inputs.subscribe{ it ->
      log.error("The expected input files or Space Ranger results files for ${it.run_id} are missing.")
    }
  

    // create tuple of [metadata, index, probeset file, fastq_dir, cytaimage file, image file]
    spaceranger_reads = spatial_channel.make_spatial
      .map{ meta_in ->
        def meta = meta_in.clone()

        // probeset logic
        def use_probeset = meta.technology in visium_probe_techs
        def probeset_file = use_probeset ? file(meta.visium_probeset) : []

        // cytaimage logic
        def cytaimage_file = getImageFiles("${meta.files_directory}/cytaimage", true)
        if (meta.technology in non_cytassist_techs) {
          if (cytaimage_file) log.error("Did not expect a cytaimage file for ${meta.technology} in ${meta.files_directory}/cytaimage but found ${cytaimage_file.size()} files.")
          cytaimage_file = []
        } else { // we have a cytassist tech
          if (cytaimage_file.size() != 1) {
            // fully error if no cytaimage, since it is required
            error("Expected exactly 1 cytaimage file in ${meta.files_directory}/cytaimage but found ${cytaimage_file.size()} files.")
            cytaimage_file = []
          } else {
            // we correctly have one cytaimage file, so index it out
            cytaimage_file = cytaimage_file[0]
          }
        } 

        // look for any additional images and check that there are not too many
        def image_file = getImageFiles("${meta.files_directory}/image")
        if (image_file.size() > 1) {
          log.error("Expected no more than 1 brightfield image file for ${meta.technology} in ${meta.files_directory}/image but found ${image_file.size()} files.")
          image_file = []
        } 
        def colorizedimage_file = getImageFiles("${meta.files_directory}/colorizedimage")
        if (colorizedimage_file.size() > 1) {
          log.error("Expected no more than 1 colorized image file for ${meta.technology} in ${meta.files_directory}/colorizedimage but found ${colorizedimage_file.size()} files.")
          colorizedimage_file = []
        } 
        def darkimage_file = getImageFiles("${meta.files_directory}/darkimage")
        if (darkimage_file.size() > 6) {
          log.error("Expected between 1-6 dark image files for ${meta.technology} in ${meta.files_directory}/darkimage but found ${darkimage_file.size()} files.")
          darkimage_file = []
        } 

        // Determine which image file to pass in to spaceranger 
        // Pair files and their types, and subset to only existing files
        def present_images = [
            [image_file, "image"],
            [colorizedimage_file, "colorizedimage"],
            [darkimage_file, "darkimage"]
          ]
          .findAll{ it[0].size() > 0 } 

        // confirm we have the right number of images:
        // no more than 1 image for any technology, but require 1 for non-cytassist
        if (present_images.size() > 1) {
          // ensure informative print statement with info about the detected files
          error("Expected no more than 1 type of provided image in ${meta.files_directory} but found: ${present_images.collect { "${it[1]}: ${it[0]}" }.join(", ")}")
        } else if (meta.technology in non_cytassist_techs && present_images.size() == 0) {
          error("At least one image file is required for ${meta.technology} but none were found in ${meta.files_directory}.")
        }

        // Simultaneously define the the image file itself and image_type (the spaceranger flag)
        def (selected_image_file, image_type) = present_images[0] ?: [[], ""]
        meta.visium_image_type = image_type

        // prepared input to spaceranger process
        [
          meta,
          file(meta.cellranger_index, type: 'dir'),
          probeset_file,
          file("${meta.files_directory}/fastq", type: 'dir'),
          cytaimage_file, 
          selected_image_file
        ]
    }
    .branch {
      hd: it[0].technology.contains("visium-hd")
      non_hd: true
    }

    // run spaceranger and mix back up
    spaceranger(spaceranger_reads.non_hd)
    spaceranger_hd(spaceranger_reads.hd)

    // gather spaceranger output for completed libraries
    // make a list of metadata (read from prior output) and prior results directory
    spaceranger_quants_ch = spatial_channel.has_spatial
      .map{ meta ->
        [
          readMeta(file("${meta.spaceranger_results_dir}/scpca-meta.json")),
          file(meta.spaceranger_results_dir, type: 'dir')
        ]
      }

    // log error about segmented_outputs as needed
    spaceranger_hd.out
      .subscribe{ meta, out_dir ->
        if (meta.visium_image_type == "image" && !file("${out_dir}/segmented_exists.txt").exists()) { 
          log.error("Expected segmented_outputs directory for ${meta.library_id} with brightfield image, but Space Ranger did not create it.")
        }
      }

    spaceranger_output_ch = spaceranger.out
      .mix(spaceranger_hd.out)
      .mix(spaceranger_quants_ch)


// generate metadata.json
spaceranger_publish(spaceranger_output_ch)

  // tuple of metadata, path to spaceranger output directory, and path to metadata json file
  emit: spaceranger_publish.out

}
