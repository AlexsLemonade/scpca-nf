
process mpileup{
  container params.BCFTOOLS_CONTAINER
  label 'cpus_2'
  tag "${meta.multiplex_run_id}"
  input:
    tuple val(meta), path(bamfiles), path(bamfiles_index)
    tuple path(ref_fasta), path(ref_index)
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
}

