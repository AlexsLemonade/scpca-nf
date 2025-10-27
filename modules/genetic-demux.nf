
// include processes
include { star_bulk } from './bulk-star.nf'
include { pileup_multibulk } from './bulk-pileup.nf'
include { starsolo_map } from './starsolo.nf'
include { cellsnp_vireo } from './cellsnp.nf'


workflow genetic_demux_vireo {
  take:
    multiplex_run_ch
    unfiltered_runs_ch
    cell_barcodes // map of cell barcode files for each technology
    bulk_techs // list of bulk technologies
  main:
    // add vireo publish directory, vireo directory, and barcode file to meta
    multiplex_ch = multiplex_run_ch
      .map{ meta_in ->
        def meta = meta_in.clone();
        meta.vireo_publish_dir = "${params.checkpoints_dir}/vireo";
        meta.vireo_dir = "${meta.vireo_publish_dir}/${meta.library_id}-vireo";
        meta.barcode_file = "${params.barcode_dir}/${cell_barcodes[meta.technology]}";
        meta // return modified meta object
      }
       // split based in whether repeat_genetic_demux is true and a previous dir exists
      .branch{ it ->
        def stored_ref_assembly = Utils.getMetaVal(file("${it.vireo_dir}/scpca-meta.json"), "ref_assembly")
        make_demux: (
          // input files exist
          it.files_directory && file(it.files_directory, type: "dir").exists() && (
            params.repeat_genetic_demux
            || !file(it.vireo_dir).exists()
            || it.ref_assembly != stored_ref_assembly
          )
        )
        has_demux: file(it.vireo_dir).exists()
        missing_inputs: true
      }

    // send run ids in multiplex_ch.missing_inputs to log
    multiplex_ch.missing_inputs
      .subscribe{ it ->
        log.error("The expected input data or vireo results files for ${it.run_id} are missing.")
      }

    // get the bulk samples that correspond to multiplexed samples
    bulk_samples = multiplex_ch.make_demux
      .map{ it -> [it.sample_id.tokenize(",")] } // split out sample ids into a tuple
      .transpose() // one element per sample id
      .collect()

    // make a channel of the bulk samples we need to process
    bulk_ch = unfiltered_runs_ch
      .filter{ it.technology in bulk_techs && it.sample_id in bulk_samples.getVal() }

    // map bulk samples
    star_bulk(bulk_ch)

    // pileup bulk samples by multiplex groups
    pileup_multibulk(multiplex_ch.make_demux, star_bulk.out)

    // map multiplexed single cell samples
    starsolo_map(multiplex_ch.make_demux, cell_barcodes)

    // call cell snps and genotype cells
    cellsnp_vireo(starsolo_map.out.bam,  starsolo_map.out.quant, pileup_multibulk.out)

    // construct demux output for skipped as [meta, vireo_dir] & join newly processed libraries
    demux_out = multiplex_ch.has_demux
      .map{ meta ->
        [
          Utils.readMeta(file("${meta.vireo_dir}/scpca-meta.json")),
          file(meta.vireo_dir, type: 'dir', checkIfExists: true)
        ]
      }
      .mix(cellsnp_vireo.out)

  emit:
    demux_out
}
