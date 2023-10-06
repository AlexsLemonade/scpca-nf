
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
    // add vireo publish directory, vireo directory, and barcode file to meta
    multiplex_ch = multiplex_run_ch
      .map{
        def meta = it.clone();
        meta.vireo_publish_dir = "${params.checkpoints_dir}/vireo";
        meta.vireo_dir = "${it.vireo_publish_dir}/${it.library_id}-vireo";
        meta.barcode_file = "${params.barcode_dir}/${params.cell_barcodes[it.technology]}";
        meta // return modified meta object
      }
       // split based in whether repeat_mapping is false and a previous dir exists
      .branch{
          has_demux: (!params.repeat_genetic_demux
                      && file(it.vireo_dir).exists()
                      && Utils.getMetaVal(file("${it.vireo_dir}/scpca-meta.json"), "ref_assembly") == "${it.ref_assembly}"
                     )
          make_demux: true
       }


    // get the bulk samples that correspond to multiplexed samples
    bulk_samples = multiplex_ch.make_demux
      .map{[it.sample_id.tokenize(",")]} // split out sample ids into a tuple
      .transpose() // one element per sample id
      .collect()

    // make a channel of the bulk samples we need to process
    bulk_ch = unfiltered_runs_ch
      .filter{it.technology in params.bulk_techs}
      .filter{it.sample_id in bulk_samples.getVal()}

    // map bulk samples
    star_bulk(bulk_ch)

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
                         file(meta.vireo_dir, type: 'dir')
                        )}
      .mix(cellsnp_vireo.out)

  emit:
    demux_out
}





