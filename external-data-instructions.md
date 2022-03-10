<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**  *generated with [DocToc](https://github.com/thlorenz/doctoc)*

- [How to use `scpca-nf` as an external user](#how-to-use-scpca-nf-as-an-external-user)
  - [Prepare the metadata file](#prepare-the-metadata-file)
  - [Configuring `scpca-nf` for your environment](#configuring-scpca-nf-for-your-environment)
    - [Configuration files](#configuration-files)
    - [Setting up a profile in the configuration file](#setting-up-a-profile-in-the-configuration-file)
    - [Using `scpca-nf` with AWS](#using-scpca-nf-with-aws)
  - [Adjust optional parameters](#adjust-optional-parameters)
  - [Special considerations for using `scpca-nf` with spatial transcriptomics libraries](#special-considerations-for-using-scpca-nf-with-spatial-transcriptomics-libraries)
  - [Output files](#output-files)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

# How to use `scpca-nf` as an external user 

In order to use `scpca-nf` to process your own data, you will need to make sure you have the following installed: 

- [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html#installation)
- [Docker](https://www.nextflow.io/docs/latest/docker.html) or [Singlularity](https://www.nextflow.io/docs/latest/singularity.html) (installed either locally or on an HPC)

You will also need to have files organized so that all the sequencing files for each run are in their own directory or folder.
Any sequencing runs that contain multiple samples must be demultiplexed and FASTQ files must be placed into separate distinct folders. 
If a run has been re-sequenced for any reason (e.g. to increase coverage), all FASTQ files should be combined into the same folder. 

Finally, you will need to create a metadata file and a nextflow configuration file (see below).
Once you have set up your environment and created these files you will be able to start your run as follows, adding any additional optional parameters that you may choose: 

```bash
nextflow run AlexsLemonade/scpca-nf \
  -r v0.2.4 \
  -config my_config.config \
  --run_metafile <path/to/metadata_file>
```

This will pull the `scpca-nf` workflow directly from Github, using the `v0.2.4` version, and run it based on the settings in the local configuration file `my_config.config`.

**Note:** `scpca-nf` is under active development.
We strongly encourage you to use a release tagged version of the workflow, set here with the `-r` flag.
Released versions can be found on the [`scpca-nf` repo releases page](https://github.com/AlexsLemonade/scpca-nf/releases).

## Prepare the metadata file 

Using `scpca-nf` requires a metadata file as a TSV (tab separated values) file, where each sequencing run to be processed is a row and columns contain associated information about that run. 

For each run, you will need to provide a Run ID (`scpca_run_id`), library ID (`scpca_library_id`), and sample ID (`scpca_sample_id`). 
Each row should contain a unique run ID, corresponding to a sequencing run or set of FASTQ files. 
If the same sequencing library was sequenced across multiple flowcells (to increase coverage, for example), these should be combined into a single row that includes all FASTQ files. 
Each run should also have a unique directory where all FASTQ files associated with that run can be found.
For example, row 1 might contain all information about a single-cell RNA-seq library and the corresponding FASTQ files. 
If that library has a corresponding CITE-seq library, and therefore a separate set of FASTQ files, the CITE-seq library should have it's own row and its own unique run ID. 

The library ID will be unique for each set of cells that have been isolated from a sample and have undergone droplet generation. 
For single-cell/single-nuclei RNA-seq runs, the library ID should be unique for each sequencing run.
For libraries that have corresponding CITE-seq, they should share the same library ID as the associated single-cell/single-nuclei RNA-seq run, indicating that the sequencing data has been generated from the same group of cells. 

Finally, the sample ID will indicate the unique tissue or source from which a sample was collected. 
If you have two libraries that have been generated from the same original tissue, then they will share the same sample ID. 

For more information on understanding the difference between library and sample IDs, see the [FAQ on library and sample IDs in the ScPCA portal documentation](https://scpca.readthedocs.io/en/latest/faq.html#what-is-the-difference-between-samples-and-libraries). 

Before using the workflow with data that you might plan to submit to ScPCA, please be sure to obtain a list of sample identifiers to use for your samples from the Data Lab. 
We will provide IDs that can be used for `scpca_run_id`, `scpca_library_id`, and `scpca_sample_id` based on the number and types of samples that are being processed to avoid overlap with existing sample identifiers. 

To run the workflow, you will need to create a tab separated values (TSV) metadata file with the following required columns: 

| column_id       | contents                                                       |
|-----------------|----------------------------------------------------------------|
| `scpca_run_id`    | A unique run ID                                              |
| `scpca_library_id`| A unique library ID for each unique set of cells             |
| `scpca_sample_id` | A unique sample ID for each tissue or unique source          |
| `technology`      | Sequencing/library technology used <br> For single-cell/single-nuclei libraries use either `10Xv2`, `10Xv2_5prime`, `10Xv3`, or `10Xv31`. <br> For CITE-seq libraries use either `CITEseq_10Xv2`, `CITEseq_10Xv3`, or `CITEseq_10Xv3.1` <br> For bulk RNA-seq use either `single_end` or `paired_end`. <br> For spatial transcriptomics use either `visium_v1` or `visium_v2`      |
| `seq_unit`        | Sequencing unit (one of: `cell`, `nucleus`, `bulk`, or `spot`)|
| `files_directory` | path/uri to directory containing fastq files (unique per run) |

The following columns may be necessary for running other data modalities (CITE-seq, spatial trancriptomics) or are optional and can be included in the metadata file if desired: 

| column_id       | contents                                                       |
|-----------------|----------------------------------------------------------------|
| `submitter_id`    | Original sample identifier defined by user (for reference only; optional)|
| `submitter`       | Name of user submitting name/id  (optional)                  |
| `feature_barcode_file`| path/uri to directory containing the feature barcode sequences (only required for CITE-seq)  |	
| `feature_barcode_geom`| A salmon `--read-geometry` layout string/ See https://github.com/COMBINE-lab/salmon/releases for details (only required for CITE-seq) |
| `slide_section`   | The slide section for spatial transcriptomics samples (only required for spatial transcriptomics) |
| `slide_serial_number`| The slide serial number for spatial transcriptomics samples (only required for spatial transcriptomics)   |
| `files`           | All sequencing files in the run folder, `;`-separated (only required for spatial transcriptomics)  |

We have provided an example metadata file for reference that can be found in [`examples/example_metadata.tsv`](examples/example_metadata.tsv).

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

Workflow parameters can also be set in a [configuration file](https://www.nextflow.io/docs/latest/config.html#configuration-file) by setting the values `params.run_metafile` and `params.outdir` as follows.

We could first create a file `my_config.config` (or a filename of your choice) with the following contents:  

```groovy
// my_config.config
params.run_metafile = '<path/to/metadata_file>'
params.outdir = '<path/to/output>'
```

This file is then used with the `-config` (or `-c`) argument at the command line:

```sh
nextflow run AlexsLemonade/scpca-nf \
  -r v0.2.4 \
  -config my_config.config 
```

For reference, we provide an example template configuration file, [user_template.config](https://github.com/AlexsLemonade/scpca-nf/blob/main/examples/user_template.config), which includes some other workflow parameters that may be useful, as well as an example of profile configuration, discussed below. 

See the [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) and the below sections for more detail on creating your own configuration file.

### Setting up a profile in the configuration file

Local running may be sufficient for small jobs or testing, but you will most likely want to run your workflow in a high performance computing environment (HPC), such as an institutional computing cluster or on a cloud service like AWS.
To do this, we recommend using [Nextflow profiles](https://www.nextflow.io/docs/latest/config.html#config-profiles) to encapsulate settings like the [`executor`](https://www.nextflow.io/docs/latest/executor.html) that will be used to run each process and associated details that may be required, such as queue names or the container engine (i.e. [Docker](https://www.nextflow.io/docs/latest/docker.html) or [Singularity](https://www.nextflow.io/docs/latest/singularity.html)) your system uses.
You will likely want to consult your HPC documentation and/or support staff to determine recommended settings.

In our example template file [user_template.config](https://github.com/AlexsLemonade/scpca-nf/blob/main/examples/user_template.config), we define a profile named `cluster` which could be invoked with the following command:

```sh
nextflow run AlexsLemonade/scpca-nf \
  -r v0.2.4 \
  -config user_template.config \
  -profile cluster
```

### Using `scpca-nf` with AWS

At the Data Lab, we use Nextflow with the [Amazon Web Services (AWS) Batch](https://aws.amazon.com/batch/) compute environment. 
If you are interested in using AWS batch with `scpca-nf`, we provide some basic instructions here to get you started.
Be aware, AWS management can be quite complex, with many interacting parts and mysterious acronyms.
We encourage you to read the official Nextflow instructions for [running pipelines on AWS](https://www.nextflow.io/docs/latest/awscloud.html), which includes information about security and permission settings that are beyond the scope of this document.

To run `scpca-nf`, you will need to set up at least one batch queue and an associated compute environment configured with a custom Amazon Machine Image (AMI) prepared according to the [Nextflow instructions](https://www.nextflow.io/docs/latest/awscloud.html#custom-ami). 
You will also need an [S3 bucket](https://aws.amazon.com/s3/) path to use as the Nextflow `work` directory for intermediate files.
As the intermediate files can get quite large, you will likely want to set up a [lifecycle rule](https://docs.aws.amazon.com/AmazonS3/latest/userguide/object-lifecycle-mgmt.html) to delete files from this location after a fixed period of time (e.g., 30 days).


In most Batch queue setups, each AWS compute node has a fixed amount of disk space. 
We found it useful to have two queues: one for general use and one for jobs that may require larger amounts of disk space. 
The two compute environments use the same AMI, but use [Launch Templates](https://docs.aws.amazon.com/batch/latest/userguide/launch-templates.html) to configure the nodes on launch with different amounts of disk space. 
Currently, our default queue is configured with a disk size of 128 GB for each node, and our "bigdisk" queue has 1000 GB of disk space.
The queue used by each process is determined by Nextflow labels and associated profile settings.

The Data Lab's [AWS Batch config file](https://github.com/AlexsLemonade/scpca-nf/blob/main/config/profile_awsbatch.config) may be helpful as a reference for creating a profile for use with AWS, but note that the queues and file locations listed there are not publicly available, so these will need to be set to different values your own profile.

## Special considerations for using `scpca-nf` with spatial transcriptomics libraries 

To process spatial transcriptomic libraries, all FASTQ files for each sequencing run and the associated `.jpg` file must be inside its own directory. 
The `run_metafile` must also contain columns with the `slide_section`, `slide_serial_number`, and list of all `files` inside the run directory, including the `.jpg` image. 

You will also need to provide a [docker image](https://docs.docker.com/get-started/) that contains the [Space Ranger software from 10X Genomics](https://support.10xgenomics.com/spatial-gene-expression/software/downloads/latest). 

After building the docker image, add the path/uri to the image to `params.SPACERANGER_CONTAINER` in the `user_template.config` file. 

## Output files 

Upon completion of the `scpca-nf` workflow, the results will be published to the specified `outdir`. 
Within the `outdir`, two folders will be present, `publish` and `internal`. 

The `publish` folder will contain the final output files produced by the workflow and the files that are typically available for download on the ScPCA portal. 

Within the `publish` folder, all files pertaining to a specific sample will be nested within a folder labeled with the sample ID.
All files in that folder will be prefixed by library ID, with the following suffixes:  `_unfiltered.rds`, `_filtered.rds`, `_metadata.json`, and `_qc.html`. 
The `_unfiltered.rds` and `_filtered.rds` files contain the quantified gene expression data as a [`SingleCellExperiment` object](https://bioconductor.org/packages/release/bioc/html/SingleCellExperiment.html).
For more information on the contents of these files, see the [ScPCA portal docs section on single cell gene expression file contents](https://scpca.readthedocs.io/en/latest/sce_file_contents.html).

See below for the expected structure of the `publish` folder: 

```
publish
└── SCPCS999991
    ├── SCPCL999991_filtered.rds
    ├── SCPCL999991_metadata.json
    ├── SCPCL999991_qc.html
    └── SCPCL999991_unfiltered.rds
```

If bulk libraries were processed, a `bulk_quant.tsv` and `bulk_metadata.tsv` summarizing the counts data and metadata across all libraries will also be present in the `publish` directory. 

The `internal` folder will contain intermediate files that are produced by individual steps of the workflow, including mapping with `salmon` and quantification with `alevin-fry`. 
The `rad` folder contains the output from running [`salmon alevin`](https://salmon.readthedocs.io/en/latest/alevin.html#alevin) with the `--rad` flag, while the `af` folder contains the outputs from [`alevin-fry`](https://alevin-fry.readthedocs.io/en/latest/index.html). 
These are created in part to allow restarting the workflow from internal checkpoints (in particular so the initial `salmon` mapping does not need to be repeated).
These also may contain log files and other outputs that can be useful for troubleshooting or alternative analyses.

All files pertaining to a specific library will be nested within a folder labeled with the library ID.
Additionally, for each run, all files related to that run will be inside a folder labeled with the run ID followed by the type of run (i.e. `rna` or `features` for CITE-seq) and nested within the library ID folder.

If bulk libraries are processed, there will be an additional `salmon` folder that contains the output from running [`salmon`](https://salmon.readthedocs.io/en/latest/file_formats.html#fileformats) on each library processed. 

## Repeat mapping 

By default, `scpca-nf` is setup to skip the mapping steps for any libraries in which the output files from the mapping step exists (i.e. the `.rad` files from `salmon alevin` and `quant.sf` files from `salmon`). 

If you would like to repeat mapping, regardless of if the files exist, you can do so by using the `--repeat_mapping` flag at the command line: 
```sh
nextflow run AlexsLemonade/scpca-nf \
  -r v0.2.4 \
  --repeat_mapping
```