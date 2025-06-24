#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process cellranger_flex_single {
  container params.CELLRANGER_CONTAINER
  publishDir "${meta.cellranger_multi_publish_dir}", mode: 'copy'
  tag "${meta.run_id}-cellranger-multi"
  label 'cpus_12'
  label 'mem_24'
  input:
    tuple val(meta), path(fastq_dir), path(cellranger_index), path(flex_probeset)
  output:
    tuple val(meta), path("${out_id}/outs/per_sample_outs/*", type: 'dir'), emit: "per_sample"
    tuple val(meta), path("${out_id}/outs/multi"), emit: "multi"
    tuple val(meta), path("${out_id}/outs/config.csv"), emit: "config"
    tuple val(meta), path("${out_id}/_versions"), emit: "versions"
  script:
    out_id = file(meta.cellranger_multi_results_dir).name
    meta += Utils.getVersions(workflow, nextflow)
    meta_json = Utils.makeJson(meta)
    """
  
    # create config file
    create_cellranger_config.py \
      --config flex_config.csv \
      --transcriptome_reference ${cellranger_index} \
      --probe_set_reference ${flex_probeset} \
      --fastq_dir ${fastq_dir}

    # run cellranger multi
    cellranger multi \
      --id=${out_id} \
      --csv=flex_config.csv \
      --localcores=${task.cpus} \
      --localmem=${task.memory.toGiga()}

    # write metadata and add cellranger version 
    echo '${meta_json}' > ${out_id}/scpca-meta.json

    # remove Cell Ranger intermediates directory
    rm -rf ${out_id}/SC_MULTI_CS
    """
  stub:
    out_id = file(meta.cellranger_multi_results_dir).name
    meta += Utils.getVersions(workflow, nextflow)
    meta_json = Utils.makeJson(meta)
    """
    mkdir -p ${out_id}/outs/per_sample_outs/${meta.sample_id}
    mkdir -p ${out_id}/outs/multi
    touch ${out_id}/outs/config.csv
    touch ${out_id}/_versions
    echo '${meta_json}' > ${out_id}/scpca-meta.json
    """
}

process cellranger_flex_multi {
  container params.CELLRANGER_CONTAINER
  publishDir "${meta.cellranger_multi_publish_dir}", mode: 'copy'
  tag "${meta.run_id}-cellranger-multi"
  label 'cpus_12'
  label 'mem_24'
  input:
    tuple val(meta), path(fastq_dir), path(cellranger_index), path(flex_probeset)
    path multiplex_pools_file
  output:
    tuple val(meta), 
      path("${out_id}/outs/per_sample_outs/*", type: 'dir'), 
      path("${out_id}/_versions"),
      path("${out_id}/outs/config.csv"),
      emit: "per_sample"
    tuple val(meta),
      path("${out_id}/outs/multi", type = 'dir'), 
      path("${out_id}/_versions"),
      path("${out_id}/outs/config.csv"),
      emit: "multi"
    tuple val(meta), path("${out_id}/outs/config.csv"), emit: "config"
    tuple val(meta), path("${out_id}/_versions"), emit: "versions"
  script:
    out_id = file(meta.cellranger_multi_results_dir).name
    meta += Utils.getVersions(workflow, nextflow)
    meta_json = Utils.makeJson(meta)
    """
  
    # create config file
    create_cellranger_config.py \
      --config flex_config.csv \
      --transcriptome_reference ${cellranger_index} \
      --probe_set_reference ${flex_probeset} \
      --fastq_dir ${fastq_dir} \
      --multiplex_pools_file ${multiplex_pools_file} \
      --library_id ${meta.library_id}

    # print for debugging
    cat flex_config.csv

    # run cellranger multi
    cellranger multi \
      --id=${out_id} \
      --csv=flex_config.csv \
      --localcores=${task.cpus} \
      --localmem=${task.memory.toGiga()}

    # write metadata
    echo '${meta_json}' > ${out_id}/scpca-meta.json

    # remove Cell Ranger intermediates directory
    rm -rf ${out_id}/SC_MULTI_CS
    """
  stub:
    out_id = file(meta.cellranger_multi_results_dir).name
    sample_ids = meta.sample_id.tokenize(",")
    meta += Utils.getVersions(workflow, nextflow)
    meta_json = Utils.makeJson(meta)
    """
    ${sample_ids.collect { "mkdir -p ${out_id}/outs/per_sample_outs/${it}" }.join("\n")}
    mkdir -p ${out_id}/outs/multi
    touch ${out_id}/outs/config.csv
    touch ${out_id}/_versions
    echo '${meta_json}' > ${out_id}/scpca-meta.json
    """
}



workflow flex_quant{
  take: 
    flex_channel // a channel with a map of metadata for each flex library to process
    flex_probesets // map of probe set files for each technology
    pool_file // path to file with barcode IDs for each sample when using multiplexed 10x flex
  main:

    flex_channel = flex_channel
      // add sample names and output directory to metadata
      .map{
        def meta = it.clone();
        meta.cellranger_multi_publish_dir =  "${params.checkpoints_dir}/cellranger-multi/${meta.library_id}";
        meta.cellranger_multi_results_dir = "${meta.cellranger_multi_publish_dir}/${meta.run_id}-cellranger-multi";
        meta.flex_probeset = "${params.probes_dir}/${flex_probesets[meta.technology]}";
        meta // return modified meta object
      }
      .branch{
        // branch based on if cellranger results exist or repeat mapping is used 
        make_cellranger_flex: (
          // input files exist
          it.files_directory && file(it.files_directory, type: 'dir').exists() && (
            params.repeat_mapping
            || !file(it.cellranger_multi_results_dir).exists()
            || Utils.getMetaVal(file("${it.cellranger_multi_results_dir}/scpca-meta.json"), "ref_assembly") != "${it.ref_assembly}"
          )
        )
        has_cellranger_flex: file(it.cellranger_multi_results_dir).exists()
        missing_inputs: true
      }

    // send run ids in flex_channel.missing_inputs to log
    flex_channel.missing_inputs
      .subscribe{
        log.error("The expected input files or cellranger multi results files for ${it.run_id} are missing.")
      }

    // tuple of inputs for running cell ranger regardless of multi or single
    // this channel only has libraries that need to be processed through cellranger
    flex_reads = flex_channel.make_cellranger_flex
      .map{ meta -> tuple(
        meta, 
        file(meta.files_directory, type: 'dir', checkIfExists: true),
        file(meta.cellranger_index, type: 'dir', checkIfExists: true),
        file(meta.flex_probeset, type: 'file', checkIfExists: true)
      )}
      .branch{ it -> 
        single: it[0]["technology"].contains("single")
        multi: it[0]["technology"].contains("multi")
      }
    
    // run cellranger flex single
    cellranger_flex_single(flex_reads.single)

    // run cellranger multiplexed
    cellranger_flex_multi(flex_reads.multi, file(pool_file))

    // get version files from output
    // first mix multi and single versions back together and then will combine with cellranger output
    versions_ch = cellranger_flex_single.out.versions.mix(cellranger_flex_multi.out.versions)
      .map{ meta, versions_file -> tuple(
        meta.library_id,
        versions_file
      )}

    // transpose cellranger multi output to have one row per output folder
    // for multiplexed data, the directory with cellranger output is in the per_sample_outs folder
    cellranger_flex_multi_flat_ch = cellranger_flex_multi.out.per_sample
      .transpose() // [meta, out_dir, versions, config]
      .map{
        def updated_meta = it[0].clone(); // clone meta before replacing sample ID
        def out_dir = it[1]; // path to individual output dir
        def sample_id = out_dir.name; // name of individual output directory is sample id
        // check that name of out directory is in expected sample IDs 
        def expected_sample_ids = updated_meta.sample_id.split(",")
        if(!(sample_id in expected_sample_ids)) {
            log.warn("${sample_id} found in output folder from cellranger multi for ${updated_meta.library_id} but does not match expected sample ids: ${expected_sample_ids}")
        }

        // update sample ID and return new meta with path to h5 file for each sample
        updated_meta.sample_id = sample_id;
        return [updated_meta, file("${out_dir}/count/sample_raw_feature_bc_matrix", type: 'dir')]
      }

    // combine single and multi outputs
    // for singleplexed data, the raw output is in the multi folder
    cellranger_flex_ch = cellranger_flex_single.out.multi
      // get path to raw output directory for singleplexed
      .map{ meta, out_dir -> tuple(
        meta, 
        file("${out_dir}/outs/multi/count/raw_feature_bc_matrix", type: 'dir')
      )}
    .mix(cellranger_flex_multi_flat_ch)
    .map{[it[0]["library_id"]] + it } // pull out library id for joining with versions
    .join(versions_ch)
    .map{it.drop(1)} // remove library_id to get [meta, out directory, versions]

    // make sure the libraries that we are skipping processing on have the correct channel format 
    // path to the raw output is dependent on single or multiplexed so split up skipped libraries based on technology
    // transpose to have one sample ID for each row
    // use existing meta to define the output directory for each sample in the library 
    has_flex_ch = flex_channel.has_cellranger_flex
      // [tuple(sample_ids), meta]
      .map{ meta -> tuple(
        meta.sample_id.tokenize(","),
        Utils.readMeta(file("${meta.cellranger_multi_results_dir}/scpca-meta.json"))
      )}
      .transpose() // [sample id, meta]
      // replace existing sample id and define path to existing directory with raw output
      .map{ sample_id, meta ->
        def updated_meta = meta.clone();
        // path depends on whether singleplex or multiplex
        def demux_h5_file
        if(meta.technology.contains("single")){
            demux_h5_file = file("${meta.cellranger_multi_results_dir}/outs/multi/count/raw_feature_bc_matrix", type: 'dir')
          } else if(meta.technology.contains("multi")) {
            demux_h5_file = file("${meta.cellranger_multi_results_dir}/outs/per_sample_outs/${sample_id}/count/sample_raw_feature_bc_matrix", type: 'dir')
          }
        def versions_file = file("${meta.cellranger_multi_results_dir}/_versions");
        updated_meta.sample_id = sample_id;
        return [updated_meta, demux_h5_file, versions_file]
      }


    // Combine single, multi, and skipped libraries
    all_flex_ch = cellranger_flex_ch.mix(has_flex_ch)

  emit: all_flex_ch

}
