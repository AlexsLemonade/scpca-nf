
// include processes
include { star_bulk } from './bulk-star.nf'
include { pileup_multibulk } from './bulk-pileup.nf'
include { starsolo_map } from './starsolo.nf'
include { cellsnp_vireo } from './cellsnp.nf'


workflow genetic_demux_vireo{
  take:
    multiplex_run_ch
    unfiltered_runs_ch
  main:
    // get all bulk samples that have fastq data
    bulk_samples = unfiltered_runs_ch
      .filter{it.technology in params.bulk_techs}
      .filter{it.files_directory && file(it.files_directory, type: "dir").exists()}
      .collect{it.sample_id}

    // add vireo publish directory, vireo directory, and barcode file to meta
    multiplex_ch = multiplex_run_ch
      .map{
        def meta = it.clone();
        meta.vireo_publish_dir = "${params.checkpoints_dir}/vireo";
        meta.vireo_dir = "${meta.vireo_publish_dir}/${meta.library_id}-vireo";
        meta.barcode_file = "${params.barcode_dir}/${params.cell_barcodes[meta.technology]}";
        meta // return modified meta object
      }
       // split based in whether repeat_genetic_demux is true and a previous dir exists
      .branch{
        make_demux: (
          // input files exist
          it.files_directory
          && file(it.files_directory, type: "dir").exists()
          && it.sample_id.tokenize(",").every{it in bulk_samples.getVal()} // all samples have fastq data
          && (
            params.repeat_genetic_demux
            || !file(it.vireo_dir).exists()
            || Utils.getMetaVal(file("${it.vireo_dir}/scpca-meta.json"), "ref_assembly") != "${it.ref_assembly}"
          )
        )
        has_demux: file(it.vireo_dir).exists()
        missing_inputs: true
      }

    // send run ids in multiplex_ch.missing_inputs to log
    multiplex_ch.missing_inputs
      .subscribe{
        log.error("The expected input data or vireo results files for ${it.run_id} are missing.")
      }

    // get the bulk samples that correspond to multiplexed samples
    bulk_multiplex_samples = multiplex_ch.make_demux
      .map{[it.sample_id.tokenize(",")]} // split out sample ids into a tuple
      .transpose() // one element per sample id
      .collect()

    // make a channel of the bulk samples we need to process
    bulk_multiplex_ch = unfiltered_runs_ch
      .filter{it.technology in params.bulk_techs}
      .filter{it.sample_id in bulk_multiplex_samples.getVal()}

    // map bulk samples
    star_bulk(bulk_multiplex_ch)

    // pileup bulk samples by multiplex groups
    pileup_multibulk(multiplex_ch.make_demux, star_bulk.out)

    // map multiplexed single cell samples
    starsolo_map(multiplex_ch.make_demux)

    // call cell snps and genotype cells
    cellsnp_vireo(starsolo_map.out.bam,  starsolo_map.out.quant, pileup_multibulk.out)

    // construct demux output for skipped as [meta, vireo_dir] & join newly processed libraries
    demux_out = multiplex_ch.has_demux
      .map{meta -> tuple(
        Utils.readMeta(file("${meta.vireo_dir}/scpca-meta.json")),
        file(meta.vireo_dir, type: 'dir', checkIfExists: true)
      )}
      .mix(cellsnp_vireo.out)

  emit:
    demux_out
}
