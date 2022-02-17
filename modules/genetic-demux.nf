
// include processes
include { star_bulk } from './bulk-star.nf'
include { mpileup } from './bcftools.nf'
include { starsolo_map } from './starsolo.nf' 
include { cellsnp; vireo } from './cellsnp.nf'

import groovy.json.JsonSlurper
import groovy.json.JsonOutput

def read_meta(path) {
  meta = new JsonSlurper().parse(file(path))
  return(meta)
}


workflow genetic_demux{
  take: 
    multiplex_ch
    unfiltered_runs_ch
  main:
    // add vireo publish directory, vireo directory, and barcode file to meta
    multiplex_ch = multiplex_ch
      .map{it.vireo_publish_dir = "${params.outdir}/internal/vireo/${it.library_id}";
           it.vireo_dir = "${it.vireo_publish_dir}/${it.run_id}-vireo"; 
           it.barcode_file = "${params.barcode_dir}/${params.cell_barcodes[it.technology]}";
           it}
       // split based in whether repeat_mapping is false and a previous dir exists
      .branch{
          has_demux: !params.repeat_gdemux && file(it.vireo_dir).exists()    
          make_demux: true
       } 
    
    
    // get the bulk samples that correspond to multiplexed samples
    bulk_samples = multiplex_ch.make_demux
      .map{[it.sample_id.tokenize("_")]} // split out sample ids into a tuple
      .transpose() // one element per sample (meta objects repeated)
      .map{it[0]} // get sample ids
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

    // construct demux output for skipped & join others
    demux_out = multiplex_ch.has_demux
      .map{[read_meta("${it.vireo_dir}/scpca-meta.json"),
            it.vireo_dir
            ]}
      .mix(cellsnp_vireo.out)
  
  emit:
    demux_out
}


workflow pileup_multibulk{
  take:
    multiplex_ch // a channel of multiplex meta objects, with sample_ids joined by `_` in their ids
    bulk_mapped_ch // output of bulk mapping: [meta, bamfile, bamfile_index]
  
  main:
    //pull sample to front of bulk_mapped_ch
    sample_bulk_ch = bulk_mapped_ch
      .map{[it[0].sample_id] + it} 

    pileup_ch = multiplex_ch 
      .map{[it.sample_id.tokenize("_"), it]} // split out sample ids into a tuple
      .transpose() // one element per sample (meta objects repeated)
      .combine(sample_bulk_ch, by: 0) // combine by individual sample ids
      .groupTuple(by: 1) // group by the multiplex run meta object
      .map{[
        [ // create a meta object for each group of samples
          sample_ids: it[0],
          multiplex_run_id: it[1].run_id,
          multiplex_library_id: it[1].library_id,
          multiplex_sample_id: it[1].sample_id,
          n_samples: it[1].sample_id.split("_").length,
          n_bulk_mapped: it[2].length,
          bulk_run_ids: it[2].collect{it.run_id},
          bulk_sample_ids: it[2].collect{it.sample_id},
          bulk_library_ids: it[2].collect{it.library_id}
        ],
        it[3], // bamfiles
        it[4]  // bamfile indexes
      ]}

    mpileup(pileup_ch, [params.ref_fasta, params.ref_fasta_index])
  
  emit:
    mpileup.out

}

workflow cellsnp_vireo {
  take: 
    starsolo_bam_ch //channel of [meta, bamfile, bam.bai]
    starsolo_quant_ch //channel of [meta, starsolo_dir]
    mpileup_vcf_ch // channel of [meta_mpileup, vcf_file]
  main:
    mpileup_ch = mpileup_vcf_ch
      .map{[it[0].multiplex_library_id] + it} // pull out library id for combining
    star_mpileup_ch = starsolo_bam_ch.map{[it[0].library_id] + it} // add library id at start
      .combine(starsolo_quant_ch.map{[it[0].library_id] + it}, by: 0) // join starsolo outs by library_id 
      .map{[it[0], it[1], it[2], it[3], it[5]]} // remove redundant meta
      .combine(mpileup_ch, by: 0) // join starsolo and mpileup by library id
      .map{it.drop(1)} // drop library id
      //result: [meta, star_bam, star_bai, star_quant, meta_mpileup, vcf_file]

    cellsnp(star_mpileup_ch) \
      | vireo
  emit:
    vireo.out
}
