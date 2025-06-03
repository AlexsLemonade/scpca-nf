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
    tuple val(meta), path(out_id)
  script:
    out_id = file(meta.cellranger_multi_publish_dir).name
    meta += Utils.getVersions(workflow, nextflow)
    meta_json = Utils.makeJson(meta)
    """
  
    # create config file
    config_file="${meta.library_id}-config.csv"

    python create_cellranger_config.py \
      --config \$config_file \
      --transcriptome_reference ${cellranger_index} \
      --probe_set_reference ${flex_probeset} \
      --sample_id ${meta.cr_sample_id} \
      --fastq_dir ${fastq_dir}

    # run cellranger multi
    cellranger multi \
      --id=${out_id} \
      --csv=\$config_file \
      --localcores=${task.cpus} \
      --localmem=${task.memory.toGiga()}

    # write metadata
    echo '${meta_json}' > ${out_id}/scpca-meta.json
    """
  stub:
    out_id = file(meta.cellranger_multi_publish_dir).name
    meta += Utils.getVersions(workflow, nextflow)
    meta_json = Utils.makeJson(meta)
    """
    mkdir -p ${out_id}/outs
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
        // get sample id in fastq files folder that is required by cellranger
        // cellranger multi only uses a single sample ID so grab from the first fastq file 
        def fastq_file = file(meta.files_directory).list().findAll{it.contains('.fastq.gz')}[0];
        meta.cr_sample_id = (fastq_file =~ /^(.+)_S.+_L.+_[R|I].+.fastq.gz$/)[0][1];
        meta.cellranger_multi_publish_dir =  "${params.checkpoints_dir}/cellranger-multi/${meta.library_id}";
        meta.flex_probeset = "${params.probes_dir}/${flex_probesets[meta.technology]}";
        meta // return modified meta object
      }
      .branch{
        make_cellranger_flex: (
          // input files exist
          it.files_directory && file(it.files_directory, type: 'dir').exists() && (
            params.repeat_mapping
            || !file(it.cellranger_multi_publish_dir).exists()
            || Utils.getMetaVal(file("${it.cellranger_multi_publish_dir}/scpca-meta.json"), "ref_assembly") != "${it.ref_assembly}"
          )
        )
        has_cellranger_flex: file(it.cellranger_multi_publish_dir).exists()
        missing_inputs: true
      }

    // send run ids in flex_channel.missing_inputs to log
    flex_channel.missing_inputs
      .subscribe{
        log.error("The expected input files or cellranger multi results files for ${it.run_id} are missing.")
      }

    // tuple of inputs for running cell ranger regardless of multi or single
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

    // TODO: Update with handling for multiplexed channel  
    // only run if a pool file is provided and there are multiplexed libraries
    multiplex_pools_ch = Channel.fromPath(pool_file)
      .splitCsv(header: true, sep: '\t')
      .map { pool_meta -> tuple(
        [pool_meta[0], pool_meta[1]],
        pool_meta
      )
      }

    flex_multi_ch = flex_reads.multi
      .map{ meta, files_dir, index -> tuple(
        meta.library_id,
        meta.sample_id.tokenize(","),
        meta,
        files_dir,
        index
      )}
      .transpose()
      .map{ library_id, sample_id, meta, files_dir, index -> tuple(
        [library_id, sample_id],
        meta,
        files_dir, 
        index
      )}
      // tuple of [[library id, sample id], meta, files_dir, index, barcode]
      .join(multiplex_pools_ch) // join by both library and sample ID


    // TODO: Run cellranger multiplexed and then join with single channel 

    // need to join back with skipped reads before outputting
    flex_quants_ch = flex_channel.has_cellranger_flex
      .map{meta -> tuple(
        Utils.readMeta(file("${meta.cellranger_multi_publish_dir}/scpca-meta.json")),
        file(meta.cellranger_multi_publish_dir, type: 'dir')
      )}

    // TODO: Make sure both single and multi are incldued before mixing
    all_flex_ch = cellranger_flex_single.out.mix(flex_quants_ch)

  emit: all_flex_ch

}
