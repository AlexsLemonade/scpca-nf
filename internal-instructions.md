# scpca-nf Data Lab Instructions

## Running scpca-nf as a Data Lab staff member

The instructions below assume that you are a member of the Data Lab with access to AWS.
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
nextflow run AlexsLemonade/scpca-nf -r v0.4.1 -profile ccdl,batch --project SCPCP000000
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
