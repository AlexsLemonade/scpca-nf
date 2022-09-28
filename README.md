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

To process single-cell and single-nuclei samples using `scpca-nf` you will need access to 24 GB of RAM and 12 CPUs, so we recommend using AWS batch.
The first step in running the workflow is ensuring that your AWS credentials are configured.

You can then run the workflow with the `batch` profile, which has been named in the `nextflow.config` file with full settings set up in the `profile_awsbatch.config` file.
Note that you will also need to use the `ccdl` profile to set file paths, and you can specify both profiles by separating them with a comma.

You can then run the workflow use the following command:

```
nextflow run AlexsLemonade/scpca-nf -profile ccdl,batch
```

When running the workflow for a project or group of samples that is ready to be released on ScPCA portal, please use the tag for the latest release:

```
nextflow run AlexsLemonade/scpca-nf -r v0.3.4 -profile ccdl,batch --project SCPCP000000
```

### Processing example data

We provide an [example of the expected outputs](./examples/README.md#example-output) after running `scpca-nf` available for external users.
If there have been major updates to the directory structure or the contents of the output, the example data should be re-processed such that the example output we provide mimics the current expected output from `scpca-nf`.

The following command should be used to run the workflow and process the example data:

```sh
nextflow run AlexsLemonade/scpca-nf \
  -profile ccdl,batch \
  --run_ids run01 \
  --run_metafile s3://scpca-references/example-data/example_metadata.tsv \
  --outdir s3://scpca-references/example-data/scpca_out
```

After successful completion of the run, the `scpca_out` folder containing the outputs from `scpca-nf` should be zipped up and stored at the following location: `s3://scpca-references/example-data/scpca_out.zip`.
Make sure to adjust the settings to make the zip file publicly accessible.
