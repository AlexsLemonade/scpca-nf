# scpca-nf

This repository holds a [Nextflow](https://www.nextflow.io) workflow to process samples as part of the ScPCA project.

Fastq files for single-cell and single-nuclei RNA-seq samples are processed using [alevin-fry](https://alevin-fry.readthedocs.io/en/latest/).
All samples are aligned, using selective alignment, to an index with transcripts corresponding to spliced cDNA and to intronic regions, denoted by alevin-fry as `splici`. 
`scpca-nf` can also process CITE-seq, bulk RNA-seq, and spatial transcriptomics. 
For more information on the processing of other modalities, please see the [ScPCA Portal docs](https://scpca.readthedocs.io/en/latest/). 


## Running scpca-nf for the ScPCA portal 

To run `scpca-nf` in its default configuration: 

```
nextflow run AlexsLemonade/scpca-nf 
```

Although running workflows locally can be done, we recommend using AWS batch for this workflow. 
The first step in running the workflow is ensuring that your AWS credentials are configured. 

You can then run the same workflow with the `batch` profile, which has been configured in the `nextflow.config` file. 

```
nextflow run AlexsLemonade/scpca-nf -profile batch
```

When running the workflow for a project or group of samples that is ready to be released on ScPCA portal, please use the tag for the latest release: 

```
nextflow run AlexsLemonade/scpca-nf -r v0.2.5 -profile batch --project SCPCP000000
```

Note that the default settings of the workflow is configured for the ALSF Childhood Cancer Data Lab computational infrastructure. 
To process samples that are not part of the ScPCA project, please see the [instructions on using `scpca-nf` with external data](external-data-instructions.md).
