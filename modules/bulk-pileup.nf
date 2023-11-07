
process mpileup{
  container params.BCFTOOLS_CONTAINER
  label 'cpus_2'
  tag "${meta.multiplex_run_id}"
  input:
    tuple val(meta), path(bamfiles), path(bamfiles_index), path(ref_fasta), path(ref_index)
  output:
    tuple val(meta), path(mpileup_file)
  script:
    mpileup_file = "${meta.multiplex_library_id}.vcf.gz"
    """
    # create sample file to use for header
    echo "${meta.sample_ids.join('\n')}" > samples.txt

    # call genotypes against reference, filter sites with missing genotypes
    # & use sample names for header (replacing file names)
    bcftools mpileup -Ou \
      --fasta-ref ${ref_fasta} \
      ${bamfiles} \
    | bcftools call -Ou --multiallelic-caller --variants-only \
    | bcftools view -Oz --genotype ^miss \
    | bcftools reheader -s samples.txt \
    > ${mpileup_file}
    """
  stub:
    mpileup_file = "${meta.multiplex_library_id}.vcf.gz"
    """
    touch ${mpileup_file}
    """
}

workflow pileup_multibulk{
  take:
    multiplex_ch // a channel of multiplex meta objects, with sample_ids joined by `,` in their ids
    bulk_mapped_ch // output of bulk mapping: [meta, bamfile, bamfile_index]

  main:
    //pull sample to front of bulk_mapped_ch
    sample_bulk_ch = bulk_mapped_ch
      .map{[it[0].sample_id] + it}

    pileup_ch = multiplex_ch
      .map{[it.sample_id.tokenize(","), it.library_id, it]} // split out sample ids into a tuple, add library_id separately
      .transpose() // one element per sample (library & meta objects repeated)
      .combine(sample_bulk_ch, by: 0) // combine by individual sample ids (use combine because of repeated values)
      .groupTuple(by: 1) // group by library id
      .map{[
        [ // create a meta object for each group of samples
          sample_ids: it[0],
          multiplex_run_id: it[2][0].run_id, // multiplex meta objects are repeated, but identical: use first element
          multiplex_library_id: it[2][0].library_id,
          multiplex_sample_id: it[2][0].sample_id,
          n_samples: it[2][0].sample_id.split(",").length,
          n_bulk_mapped: it[3].length,
          bulk_run_ids: it[3].collect{it.run_id},
          bulk_sample_ids: it[3].collect{it.sample_id},
          bulk_library_ids: it[3].collect{it.library_id}
        ],
        it[4], // bamfiles
        it[5], // bamfile indexes
        file(it[2][0].ref_fasta),
        file(it[2][0].ref_fasta_index)
      ]}

    mpileup(pileup_ch)

  emit:
    mpileup.out

}
