
// include processes
include { star_bulk } from './bulk-star.nf'
include { starsolo_map } from './starsolo.nf' 
include { pileup_multibulk } from './sambcftools.nf'
include { cellsnp_vireo } from './cellsnp.nf'


workflow genetic_demux{
  take: 
    multiplex_ch
    unfiltered_runs_ch
  main:
    // get the bulk samples that correspond to multiplexed samples
    bulk_samples = multiplex_ch
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
    pileup_multibulk(multiplex_ch, star_bulk.out)

    // map multiplexed single cell samples
    starsolo_map(multiplex_ch)

    // call cell snps and genotype cells 
    cellsnp_vireo(starsolo_map.out.bam,  starsolo_map.out.quant, pileup_multibulk.out)
  
  emit:
    cellsnp_vireo.out
}
