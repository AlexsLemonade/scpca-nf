<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**  *generated with [DocToc](https://github.com/thlorenz/doctoc)*

- [Prepare the metadata file](#prepare-the-metadata-file)
- [Configuring `scpca-nf` for your environment](#configuring-scpca-nf-for-your-environment)
  - [Configuration files](#configuration-files)
  - [Setting up a profile in the configuration file](#setting-up-a-profile-in-the-configuration-file)
  - [Using `scpca-nf` with AWS](#using-scpca-nf-with-aws)
- [Adjust optional parameters](#adjust-optional-parameters)
- [Special considerations for using `scpca-nf` with spatial transcriptomics libraries](#special-considerations-for-using-scpca-nf-with-spatial-transcriptomics-libraries)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

# How to use `scpca-nf` as an external user 

In order to use `scpca-nf` to process your own data, you will need to make sure you have the following installed: 

- [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html#installation)
- [Docker](https://docs.docker.com/get-started/#download-and-install-docker) (installed either locally or on an HPC)

You will also need to create a metadata file and a Nextflow configuration file (see below).
Once you have set up your environment and created these files you will be able to start your run as follows, adding any additional optional parameters that you may choose: 

```bash
nextflow run AlexsLemonade/scpca-nf \
  -r v0.2.4 \
  -config my_config.config \
  --run_metafile <path/to/metadata_file>
```

This will pull the `scpca-nf` workflow directly from Github, using the `v0.2.4` version, and run it based on the settings in the local configuration file `my_config.config`.

**Note:** `spca-nf` is under active development.
We strongly encourage you to use a release tagged version of the workflow, set here with the `-r` flag.
Released versions can be found on the [`scpca-nf` repo releases page](https://github.com/AlexsLemonade/scpca-nf/releases).

## Prepare the metadata file 

Using `scpca-nf` requires a metadata file where each library to be processed is a row and columns contain associated information about that library. 
You will need the following columns to be present in your metadata file... 

## Configuring `scpca-nf` for your environment

Two workflow parameters are *required* for running `scpca-nf` on your own data:

- `run_metafile`: the metadata file with sample information, prepared according to the directions above
- `outdir`: the output directory where results will be stored

By default, the workflow is set up to run in a local environment, and these parameters can be set at the command line as follows:

```sh
nextflow run AlexsLemonade/scpca-nf \
  -r v0.2.4 \
  --run_metafile <path/to/metadata_file> \
  --outdir <path/to/output>
```

Note that *workflow* parameters such as `--run_metafile` and `--outdir` are denoted at the command line with double hyphen prefix, while options that affect Nextflow itself have only a single hyphen. 

### Configuration files

Workflow parameters can also be set in a configuration file by setting the values `params.run_metafile` and `params.outdir` as follows: 

```groovy
// my_config.config
params.run_metafile = '<path/to/metadata_file>'
params.outdir = '<path/to/output>'
```

This is then used with the `-config` (or `-c`) argument at the command line:

```sh
nextflow run AlexsLemonade/scpca-nf \
  -r v0.2.4 \
  -config my_config.config 
```

We provide an example template configuration file, [user_template.config](https://github.com/AlexsLemonade/scpca-nf/blob/main/examples/user_template.config), for reference, which includes some other workflow parameters that may be useful, as well as an example of profile configuration, discussed below. 

See the [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) and the below sections for more detail on creating your own configuration file.

### Setting up a profile in the configuration file

Local running may be sufficient for small jobs or testing, but you will most likely want to run your workflow in a high performance computing environment (HPC), such as an institutional computing cluster or on a cloud service like AWS.
To do this, we recommend using [Nextflow profiles](https://www.nextflow.io/docs/latest/config.html#config-profiles) to encapsulate settings like the [`executor`](https://www.nextflow.io/docs/latest/executor.html) that will be used to run each process and associated details that may be required, such as queue names or the container engine your system uses. 

In our example template file [user_template.config](https://github.com/AlexsLemonade/scpca-nf/blob/main/examples/user_template.config), we define a profile named `cluster` which could be invoked with the following command:

```sh
nextflow run AlexsLemonade/scpca-nf \
  -r v0.2.4 \
  -config user_template.config \
  -profile cluster
```

### Using `scpca-nf` with AWS
Instructions on setting up AWS

## Adjust optional parameters

Include information on all parameters that can be altered. 

## Special considerations for using `scpca-nf` with spatial transcriptomics libraries 

Instructions on creating your own spaceranger docker image
