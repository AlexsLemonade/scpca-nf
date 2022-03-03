<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**  *generated with [DocToc](https://github.com/thlorenz/doctoc)*

- [How to use `scpca-nf` as an external user](#how-to-use-scpca-nf-as-an-external-user)
  - [Prepare the metadata file](#prepare-the-metadata-file)
  - [Configuring `scpca-nf` for your environment](#configuring-scpca-nf-for-your-environment)
    - [Setting up a profile in the config file](#setting-up-a-profile-in-the-config-file)
    - [Using `scpca-nf` with AWS Batch](#using-scpca-nf-with-aws-batch)
  - [Adjust optional parameters](#adjust-optional-parameters)
  - [Special considerations for using `scpca-nf` with spatial transcriptomics libraries](#special-considerations-for-using-scpca-nf-with-spatial-transcriptomics-libraries)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

## How to use `scpca-nf` as an external user 

In order to use `scpca-nf` to process your own data, you will need to make sure you have the following installed: 

- [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html#installation)
- [Docker](https://docs.docker.com/get-started/#download-and-install-docker) (installed either locally or on an HPC)

You will also need to create a metadata file and a Nextflow configuration file (see below).
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
The config file will allow you to define the system where you would like to execute Nextflow, allocate resources, and define parameters for the workflow. 
After creating a config file, `my_config.config`, that file can be used in conjunction with `scpca-nf` at the command line by using `-c my_config.config`. 
See the [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) and the below sections for more detail on creating your own config file. 

#### Setting up a profile in the config file
Instructions on adding a new profile


#### Using `scpca-nf` with AWS Batch

At the Data Lab, we use Nextflow with the [Amazon Web Services (AWS) Batch](https://aws.amazon.com/batch/) compute environment, so we can provide a bit more instruction related to that setup, if that is a direction you want to go. 
Be aware, AWS management can be quite complex, with many interacting parts and mysterious acronyms.
We encourage you to read the official Nextflow instructions for [running pipelines on AWS](https://www.nextflow.io/docs/latest/awscloud.html), which includes information about security and permission settings that are beyond the scope of this document.

To run `scpca-nf`, you will need to set up at least one batch queue and an associated compute environment configured with a custom Amazon Machine Image (AMI) prepared according to the [Nextflow instructions](https://www.nextflow.io/docs/latest/awscloud.html#custom-ami). 
You will also need an [S3 bucket](https://aws.amazon.com/s3/) path to use as the Nextflow `work` directory for intermediate files.
As the intermediate files can get quite large, you will likely want to set up a [lifecycle rule](https://docs.aws.amazon.com/AmazonS3/latest/userguide/object-lifecycle-mgmt.html) to delete files from this location after a fixed period of time (e.g., 30 days).

Our [AWS Batch config file](https://github.com/AlexsLemonade/scpca-nf/blob/main/config/profile_awsbatch.config) may be helpful as a starting reference, but note that the queues and file locations listed there are not publicly available, so these will need to be set in your own profile.

In most Batch queue setups, each AWS compute node has a fixed amount of disk space. 
We found it useful to have two queues: one for general use and one for jobs that may require larger amounts of disk space. 
The two compute environments use the same AMI, but use [Launch Templates](https://docs.aws.amazon.com/batch/latest/userguide/launch-templates.html) to configure the nodes on launch with different amounts of disk space. 
Currently, our default queue is configured with a disk size of 128 GB for each node, and our "bigdisk" queue has 1000 GB of disk space.
The queue used by each process is determined by Nextflow labels, in combination with settings in our [aws batch profile](https://github.com/AlexsLemonade/scpca-nf/blob/main/config/profile_awsbatch.config).


### Adjust optional parameters

Include information on all parameters that can be altered. 

### Special considerations for using `scpca-nf` with spatial transcriptomics libraries 

Instructions on creating your own spaceranger docker image
