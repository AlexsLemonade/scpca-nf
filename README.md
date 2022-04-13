# scpca-nf

This repository holds a [Nextflow](https://www.nextflow.io) workflow to process samples as part of the [Single-cell Pediatric Cancer Atlas (ScPCA) project](https://scpca.alexslemonade.org/).

Fastq files for single-cell and single-nuclei RNA-seq samples are processed using [alevin-fry](https://alevin-fry.readthedocs.io/en/latest/).
All samples are aligned, using selective alignment, to an index with transcripts corresponding to spliced cDNA and to intronic regions, denoted by alevin-fry as `splici`. 
`scpca-nf` can also process CITE-seq, bulk RNA-seq, and spatial transcriptomics. 
For more information on the processing of all modalities, please see the [ScPCA Portal docs](https://scpca.readthedocs.io/en/latest/). 

## Using scpca-nf to process your samples

The `scpca-nf` workflow is currently set up to process samples as part of [the ScPCA portal](https://scpca.alexslemonade.org/) and requires access to AWS through the Data Lab.  
For all other users, `scpca-nf` can be set up to process your samples in your computing environment by following the [instructions for using `scpca-nf` with external data](external-data-instructions.md). 

:warning: Please note that processing single-cell and single-nuclei RNA-seq samples, requires access to a high performance computing (HPC) environment that can accomodate up to 24 GB of RAM and 12 CPUs. 

To run `scpca-nf` on your own samples, you will need to complete the following steps: 

1. [Organize your files](./external-data-instructions.md#file-organization) so that each folder contains fastq files relevant to a single sequencing run. 
2. [Prepare a metadata file](./external-data-instructions.md#prepare-the-metadata-file) with one row per library containing all information needed to process your samples.
3. Set up a [configuration file](./external-data-instructions.md#configuration-files), including the [definition of a profile](./external-data-instructions.md#setting-up-a-profile-in-the-configuration-file), dictating where nextflow should execute the workflow. 

You may also [test your configuration file using example data](examples/README.md).

## Running scpca-nf as a Data Lab staff member

The instructions below assume that you are a member of the CCDL with access to AWS.
Most of the workflow settings described are configured for the ALSF Childhood Cancer Data Lab computational infrastructure. 
To process samples that are not part of the ScPCA project, please see the [instructions on using `scpca-nf` with external data](external-data-instructions.md). 

To run `scpca-nf` in its default configuration for the CCDL, you can use the following command: 

```
nextflow run AlexsLemonade/scpca-nf -profile ccdl
```
Although running workflows locally can be done, we recommend using AWS batch for this workflow. 
The first step in running the workflow is ensuring that your AWS credentials are configured. 

You can then run the same workflow with the `batch` profile, which has been configured in the `nextflow.config` file. 
Note that you will still need the `ccdl` profile, and you can specify both with by separating them with a comma. 

```
nextflow run AlexsLemonade/scpca-nf -profile ccdl,batch
```

When running the workflow for a project or group of samples that is ready to be released on ScPCA portal, please use the tag for the latest release: 

```
nextflow run AlexsLemonade/scpca-nf -r v0.2.7 -profile ccdl,batch --project SCPCP000000
```


