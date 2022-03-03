<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**  *generated with [DocToc](https://github.com/thlorenz/doctoc)*

- [How to use `scpca-nf` as an external user](#how-to-use-scpca-nf-as-an-external-user)
  - [Prepare the metadata file](#prepare-the-metadata-file)
  - [Configuring `scpca-nf` for your environment](#configuring-scpca-nf-for-your-environment)
    - [Setting up a profile in the config file](#setting-up-a-profile-in-the-config-file)
    - [Adding parameters to the config file](#adding-parameters-to-the-config-file)
    - [Using `scpca-nf` with AWS](#using-scpca-nf-with-aws)
  - [Adjust optional parameters](#adjust-optional-parameters)
  - [Special considerations for using `scpca-nf` with spatial transcriptomics libraries](#special-considerations-for-using-scpca-nf-with-spatial-transcriptomics-libraries)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

## How to use `scpca-nf` as an external user 

In order to use `scpca-nf` to process your own data, you will need to make sure you have the following installed: 

- [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html#installation)
- [Docker](https://docs.docker.com/get-started/#download-and-install-docker) (installed either locally or on an HPC)

You will also need to create a metadata file and a nextflow configuration file (see below).
Once you have set up your environment and created these files you will be able to start your run as follows, adding any additional optional parameters that you may choose: 

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

### Configuring `scpca-nf` for your environment

By default, the workflow is set up to run in a local environment. 
If choosing to run the workflow locally, you will only need to provide the path to the metadata file and output directory where you would like the results to be stored. 

```
nextflow run AlexsLemonade/scpca-nf \
 -r 0.2.3 \
  --run_metafile <path/to/metadata file> \
  --outdir <path/to/output>
```

However, if you need to run your workflow on an HPC, on AWS, or other workspace we encourage you to create your own config file.  
The config file will allow you to define the system where you would like to execute nextflow, allocate resources, and define parameters for the workflow. 
After creating a config file, `my_config.config`, that file can be used in conjunction with `scpca-nf` at the command line by using `-c my_config.config`. 
See the [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) and the below sections for more detail on creating your own config file. 

#### Setting up a profile in the config file
Instructions on adding a new profile

#### Adding parameters to the config file
Adding parameters to the config file. 

#### Using `scpca-nf` with AWS
Instructions on setting up AWS

### Adjust optional parameters

Include information on all parameters that can be altered. 

### Special considerations for using `scpca-nf` with spatial transcriptomics libraries 

Instructions on creating your own spaceranger docker image
