## scpca-nf

This repo holds the scripts to run the nextflow workflow to process samples as part of the ScPCA project. Fastq files for each sample are processed using [alevin-fry](https://alevin-fry.readthedocs.io/en/latest/). All samples are aligned, using selective alignment, to an index with transcripts corresponding to spliced cDNA and to intronic regions, denoted by alevin-fry as `splici`. 

To run this workflow use: 

```
nextflow run AlexsLemonade/scpca-nf 
```

Although running workflows locally can be done, we recommend using AWS batch for this workflow. 
The first step in running the workflow is ensuring that your AWS credentials are configured. 

You can then run the same workflow with the `batch` profile, which has been configured in the `nextflow.config` file. 

```
nextflow run AlexsLemonade/scpca-nf -profile batch
```