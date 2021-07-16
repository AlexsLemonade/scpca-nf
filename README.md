## scpca-nf

This repo holds the scripts to run the nextflow workflow to process samples as part of the ScPCA project. Fastq files for each sample are processed using [alevin-fry](https://alevin-fry.readthedocs.io/en/latest/). All samples are aligned, using selective alignment, to an index with transcripts corresponding to spliced cDNA and to intronic regions, denoted by alevin-fry as `splici`. 

To run this workflow use: 

```
nextflow run AlexsLemonade/scpca-nf 
```

This workflow can also be run using AWS batch. For instructions on how to do so see the [alsf-scpca repo.](https://github.com/AlexsLemonade/alsf-scpca#readme)