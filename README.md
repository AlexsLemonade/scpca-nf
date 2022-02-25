<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**  *generated with [DocToc](https://github.com/thlorenz/doctoc)*

- [Running scpca-nf for the ScPCA portal](#running-scpca-nf-for-the-scpca-portal)
- [How to use `scpca-nf` as an external user](#how-to-use-scpca-nf-as-an-external-user)
  - [Prepare the metadata file](#prepare-the-metadata-file)
  - [Identify where to run scpca-nf](#identify-where-to-run-scpca-nf)
    - [Using `scpca-nf` with AWS](#using-scpca-nf-with-aws)
    - [Adding a profile to `scpca-nf`](#adding-a-profile-to-scpca-nf)
  - [Choose an index](#choose-an-index)
  - [Adjust optional parameters](#adjust-optional-parameters)
  - [Special considerations for using `scpca-nf` with spatial transcriptomics libraries](#special-considerations-for-using-scpca-nf-with-spatial-transcriptomics-libraries)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

## Running scpca-nf for the ScPCA portal 

This repo holds the scripts to run the nextflow workflow to process samples as part of the ScPCA project.
Fastq files for single-cell and single-nuclei RNA-seq samples are processed using [alevin-fry](https://alevin-fry.readthedocs.io/en/latest/). All samples are aligned, using selective alignment, to an index with transcripts corresponding to spliced cDNA and to intronic regions, denoted by alevin-fry as `splici`. 
`scpca-nf` can also process CITE-seq, bulk RNA-seq, and spatial transcriptomics. 
For more information on the processing of other modalities, please see the [ScPCA Portal docs](https://scpca.readthedocs.io/en/latest/). 

To run `scpca-nf` in it's default configuration: 

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
nextflow run AlexsLemonade/scpca-nf -r v0.2.3 -profile batch --project SCPCP000000
```

## How to use `scpca-nf` as an external user 

In order to use `scpca-nf` to process your own data, you will need to follow the below steps. 

Once you have completed the below steps you will be able to setup your run as follows, adding any additional optional parameters that you may choose: 

```
nextflow run AlexsLemonade/scpca-nf \
 -r 0.2.3 \
 -c my_config.config \
  --run_metafile <path/to/metadata file> \
  --outdir <path/to/output>
```

### Prepare the metadata file 

Using `scpca-nf` requires a metadata file where each library to be processed is a row and columns contain associated information about that library. 
You will need the following columns to be present in your metadata file... 

### Identify where to run scpca-nf

The `scpca-nf` workflow is currently set up to run either on AWS or locally. 
If you would like to run using AWS, you will need to... 
However, if you need to run your workflow on an HPC or other workspace you will need to add an additional profile to the nextflow configuration files. 

#### Using `scpca-nf` with AWS
Instructions on setting up AWS

#### Adding a profile to `scpca-nf`
Instructions on adding a new profile

### Choose an index 

We have provided pre-built index files that are stored in AWS S3 storage that can be used regardless of how you choose to run the workflow.

However, if you would like to create your own index, you then must provide the directory containing the index to `--ref_dir` at the command line.  

For single-cell/single-nuclei samples you will need to provide a [`splici` index](https://combine-lab.github.io/alevin-fry-tutorials/2021/improving-txome-specificity/).
For more information on how to build a `splici` index... 

Include information about each of the files that we need and the file structure for the reference directory. 

### Adjust optional parameters

Include information on all parameters that can be altered. 

### Special considerations for using `scpca-nf` with spatial transcriptomics libraries 

Instructions on creating your own spaceranger docker image

