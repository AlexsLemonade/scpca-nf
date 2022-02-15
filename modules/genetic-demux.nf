
// include processes
include { star_bulk } from './bulk-star.nf'
include { mpileup } from './bcftools.nf'
include { starsolo_map } from './starsolo.nf' 
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
        [ // create a meta object for each group of files
          sample_ids: it[0],
          multiplex_run_id: it[1].run_id,
          multiplex_library_id: it[1].library_id,
          multiplex_sample_id: it[1].sample_id,
          n_samples: it[1].sample_id.split("_").length,
          n_bulk_mapped: it[2].length,
          bulk_run_ids: it[2].collect{it.run_id},
          bulk_run_prefixes: it[2].collect{it.files_directory}
        ],
        it[3], // bamfiles
        it[4]  // bamfile indexes
      ]}

    mpileup(pileup_ch, [params.ref_fasta, params.ref_fasta_index])
  
  emit:
    mpileup.out

}
