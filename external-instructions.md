# How to use `scpca-nf` as an external user

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->

- [Overview](#overview)
- [File organization](#file-organization)
- [Prepare the metadata file](#prepare-the-metadata-file)
- [Configuring `scpca-nf` for your environment](#configuring-scpca-nf-for-your-environment)
  - [Configuration files](#configuration-files)
  - [Setting up a profile in the configuration file](#setting-up-a-profile-in-the-configuration-file)
  - [Using `scpca-nf` with AWS Batch](#using-scpca-nf-with-aws-batch)
  - [Using `scpca-nf` on nodes without direct internet access](#using-scpca-nf-on-nodes-without-direct-internet-access)
    - [Additional reference files](#additional-reference-files)
    - [Downloading container images](#downloading-container-images)
- [Repeating mapping steps](#repeating-mapping-steps)
- [Special considerations for specific data types](#special-considerations-for-specific-data-types)
  - [Libraries with additional feature data (CITE-seq or cellhash)](#libraries-with-additional-feature-data-cite-seq-or-cellhash)
  - [Multiplexed (cellhash) libraries](#multiplexed-cellhash-libraries)
  - [Spatial transcriptomics libraries](#spatial-transcriptomics-libraries)
- [Output files](#output-files)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->


 ## Overview

Using `scpca-nf` to process your own single-cell and single-nuclei RNA-seq data requires access to a high performance computing (HPC) environment that can accomodate up to 24 GB of RAM and 12 CPUs.
After identifying the system that you will use to execute the nextflow workflow, you will need to follow the steps outlined in this document to complete the set up process.
Here we provide an overview of the steps you will need to complete:

1. **Install the necessary dependencies.**
You will need to make sure you have the following software installed on your HPC where you plan to execute the workflow:
    - [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html#installation), the main workflow engine that scpca-nf relies on.
    This can be downloaded and installed by any user, with minimal external requirements.
    - [Docker](https://docs.docker.com/get-docker/) or [Singularity](https://sylabs.io/guides/3.0/user-guide/installation.html#installation), which allows the use of container images that encapsulate other dependencies used by the workflow reproducibly.
    These usually require installation by system administrators, but most HPC systems have one available (usually Singularity).
    - Other software dependencies, as well as the workflow files themselves, are handled by Nextflow, which will download Docker or Singularity images as required.
    The `scpca-nf` workflow does not need to be downloaded separately.
    However, if nodes on your HPC do no not have direct internet access, you will need to follow [our instructions to download reference files and container images](#using-scpca-nf-on-nodes-without-direct-internet-access).

2. **Organize your files.**
You will need to have your files organized in a particular manner so that each folder contains only the FASTQ files that pertain to a single library.
See the [section below on file organization](#file-organization) for more information on how to set up your files.

3. **Create a [metadata file](#prepare-the-metadata-file).**
Create a TSV (tab-separated values) file with one sequencing library per row and pertinent information related to that sequencing run in each column.
See the [section below on preparing a metadata file](#prepare-the-metadata-file) for more information on creating a metadata file for your samples.

4. **Create a [configuration file](#configuration-files) and [define a profile](#setting-up-a-profile-in-the-configuration-file).**
Create a configuration file that stores user defined parameters and a profile indicating the system and other system related settings to use for executing the workflow.
See the [section below on configuring `scpca-nf` for your environment](#configuring-scpca-nf-for-your-environment) for more information on setting up the configuration files to run Nextflow on your system.

The standard configuration the `scpca-nf` workflow expects that compute nodes will have direct access to the internet, and will download reference files and container images with ane required software as required.
If your HPC system does not allow internet access from compute nodes, you will need to download the required reference files and software before running, [following the instructions we have provided](#using-scpca-nf-on-nodes-without-direct-internet-access).

Once you have set up your environment and created the metadata and configuration files, you will be able to start your run as follows, adding any additional optional parameters that you may choose:

```sh
nextflow run AlexsLemonade/scpca-nf \
  -config <path to config file>  \
  -profile <name of profile>
```

Where `<path to config file>` is the **relative** path to the [configuration file](#configuration-files) that you have setup and `<name of profile>` is the name of the profile that you chose when [creating a profile](#setting-up-a-profile-in-the-configuration-file).
This command will pull the `scpca-nf` workflow directly from Github, using the `v0.4.0` version, and run it based on the settings in the configuration file that you have defined.

**Note:** `scpca-nf` is under active development.
Using the above command will run the workflow from the `main` branch of the workflow repository.
To update to the latest released version you can run `nextflow pull AlexsLemonade/scpca-nf` before the `nextflow run` command.

To  be sure that you are using a consistent version, you can specify use of a release tagged version of the workflow, set below with the `-r` flag.
Released versions can be found on the [`scpca-nf` repo releases page](https://github.com/AlexsLemonade/scpca-nf/releases).

```sh
nextflow run AlexsLemonade/scpca-nf \
  -r v0.4.0 \
  -config <path to config file>  \
  -profile <name of profile>
```

For each library that is successfully processed, the workflow will return quantified gene expression data as a `SingleCellExperiment` object stored in an RDS file along with a summary HTML report and any relevant intermediate files.
For a complete description of the expected output files, see the section describing [output files](#output-files).

## File organization

You will need to have files organized so that all the sequencing files for each library are in their own directory or folder.
Each folder should be named with a unique ID, corresponding to the [`scpca_run_id` column of the metadata file](#prepare-the-metadata-file).
Any sequencing runs that contain multiple libraries must be demultiplexed and FASTQ files must be placed into separate distinct folders, with distinct run IDs as the folder name.
If the same sequencing library was sequenced across multiple flowcells (e.g., to increase coverage), all FASTQ files should be combined into the same folder.
If a library has a corresponding CITE-seq library and therefore a separate set of FASTQ files, the FASTQ files corresponding to the CITE-seq library should be in their own folder, with a unique run ID.

## Prepare the metadata file

Using `scpca-nf` requires a metadata file as a TSV (tab separated values) file, where each sequencing run to be processed is a row and columns contain associated information about that run.

For each sequencing run, you will need to provide a Run ID (`scpca_run_id`), library ID (`scpca_library_id`), and sample ID (`scpca_sample_id`).
The run ID will correspond to the name of the folder that contains the FASTQ files associated with the sequencing run.
See [the section on file organization above for more information](#file-organization).

The library ID will be unique for each set of cells that have been isolated from a sample and have undergone droplet generation.
For single-cell/single-nuclei RNA-seq runs, the library ID should be unique for each sequencing run.
For libraries that have corresponding CITE-seq or cellhash runs, they should share the same library ID as the associated single-cell/single-nuclei RNA-seq run, indicating that the sequencing data has been generated from the same group of cells.

Finally, the sample ID will indicate the unique tissue or source from which a sample was collected.
If you have two libraries that have been generated from the same original tissue, then they will share the same sample ID.

For more information on understanding the difference between library and sample IDs, see the [FAQ on library and sample IDs in the ScPCA portal documentation](https://scpca.readthedocs.io/en/latest/faq.html#what-is-the-difference-between-samples-and-libraries).

**Before using the workflow with data that you might plan to submit to ScPCA, please be sure to obtain a list of sample identifiers to use for your samples from the Data Lab.**
We will provide IDs that can be used for `scpca_run_id`, `scpca_library_id`, and `scpca_sample_id` based on the number and types of samples that are being processed to avoid overlap with existing sample identifiers.

To run the workflow, you will need to create a tab separated values (TSV) metadata file with the following required columns:

| column_id       | contents                                                       |
|-----------------|----------------------------------------------------------------|
| `scpca_run_id`    | A unique run ID                                              |
| `scpca_library_id`| A unique library ID for each unique set of cells             |
| `scpca_sample_id` | A unique sample ID for each tissue or unique source. <br> For multiplexed libraries, separate multiple samples with semicolons (`;`)          |
| `scpca_project_id` | A unique ID for each group of related samples. All results for samples with the same project ID will be returned in the same folder labeled with the project ID. |
| `technology`      | Sequencing/library technology used <br> For single-cell/single-nuclei libraries use either `10Xv2`, `10Xv2_5prime`, `10Xv3`, or `10Xv31`. <br> For CITE-seq libraries use either `CITEseq_10Xv2`, `CITEseq_10Xv3`, or `CITEseq_10Xv3.1` <br> For cellhash libraries use either `cellhash_10Xv2`, `cellhash_10Xv3`, or `cellhash_10Xv3.1` <br> For bulk RNA-seq use either `single_end` or `paired_end`. <br> For spatial transcriptomics use `visium`      |
| `seq_unit`        | Sequencing unit (one of: `cell`, `nucleus`, `bulk`, or `spot`)|
| `files_directory` | path/uri to directory containing fastq files (unique per run) |

The following columns may be necessary for running other data modalities (CITE-seq, spatial trancriptomics) or are optional and can be included in the metadata file if desired:

| column_id       | contents                                                       |
|-----------------|----------------------------------------------------------------|
| `feature_barcode_file` | path/uri to file containing the feature barcode sequences (only required for CITE-seq and cellhash samples)  |
| `feature_barcode_geom` | A salmon `--read-geometry` layout string. <br> See https://github.com/COMBINE-lab/salmon/releases/tag/v1.4.0 for details (only required for CITE-seq and cellhash samples) |
| `slide_section`   | The slide section for spatial transcriptomics samples (only required for spatial transcriptomics) |
| `slide_serial_number`| The slide serial number for spatial transcriptomics samples (only required for spatial transcriptomics)   |

We have provided an example metadata file for reference that can be found in [`examples/example_metadata.tsv`](examples/example_metadata.tsv).

## Configuring `scpca-nf` for your environment

Two workflow parameters are required for running `scpca-nf` on your own data:

- `run_metafile`: the metadata file with sample information, prepared according to the directions above.
  - This has a default value of `run_metadata.tsv`, but you will likely want to set your own file path.
- `outdir`: the output directory where results will be stored.
  - The default output is `scpca_out`, but again, you will likely want to customize this.

These parameters can be set at the command line using `--run_metafile <path to run_metafile>` or `--outdir <path to output>`, but we encourage you to set them in the configuration file, following the [configuration file setup instructions below](#configuration-files).

Note that *workflow* parameters such as `--run_metafile` and `--outdir` are denoted at the command line with double hyphen prefix, while options that affect Nextflow itself have only a single hyphen.

### Configuration files

Workflow parameters can also be set in a [configuration file](https://www.nextflow.io/docs/latest/config.html#configuration-file) by setting the values `params.run_metafile` and `params.outdir` as follows.

We could first create a file `my_config.config` (or a filename of your choice) with the following contents:

```groovy
// my_config.config
params.run_metafile = '<path to run_metafile>'
params.outdir = '<path to output>'
```

This file is then used with the `-config` (or `-c`) argument at the command line:

```sh
nextflow run AlexsLemonade/scpca-nf \
  -config my_config.config
```

For reference, we provide an example template configuration file, [`user_template.config`](examples/user_template.config), which includes some other workflow parameters that may be useful, as well as an example of configuring a profile for executing the workflow on a cluster, discussed below.

**Note:** This example tells Nextflow to use the configuration set up in the configuration file, but it does not invoke a specific profile, and will use the `standard` profile.
Under the `standard` profile, Nextflow will attempt to run the workflow locally using Docker.
This will most likely result in an error unless the minimum computing requirements (24 GB of RAM and 12 CPUs) are met on the local machine.
For more on [creating and using a profile](#setting-up-a-profile-in-the-configuration-file) see the section below.

See the [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) and the below sections for more detail on creating your own configuration file.

### Setting up a profile in the configuration file

Processing single-cell and single-nuclei samples requires access to 24 GB of RAM and 12 CPUs so you will most likely want to run your workflow in a high performance computing environment (HPC), such as an institutional computing cluster or on a cloud service like AWS.
To do this, we recommend using [Nextflow profiles](https://www.nextflow.io/docs/latest/config.html#config-profiles) to encapsulate settings like the [`executor`](https://www.nextflow.io/docs/latest/executor.html) that will be used to run each process and associated details that may be required, such as queue names or the container engine (i.e., [Docker](https://www.nextflow.io/docs/latest/docker.html) or [Singularity](https://www.nextflow.io/docs/latest/singularity.html)) your system uses.
You will likely want to consult your HPC documentation and/or support staff to determine recommended settings.

**Note:** To use the default index files, which are stored on S3, compute nodes must have access to the internet.
You may also need to supply AWS credentials for S3 access, or set `aws.client.anonymous = true` within the Nextflow profile.

In our example template file [`user_template.config`](examples/user_template.config), we define a profile named `cluster` which could be invoked with the following command:

```sh
nextflow run AlexsLemonade/scpca-nf \
  -config user_template.config \
  -profile cluster
```

### Using `scpca-nf` with AWS Batch

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

### Using `scpca-nf` on nodes without direct internet access

Some HPC systems limit the network traffic of compute nodes for security reasons.
The standard configuration of the `scpca-nf`, however, expects that reference files and container images (for docker or singularity) can be downloaded as needed.
If your system does not allow direct internet access, you will need to pre-download the required reference files to a local directory and adjust parameters to direct the workflow to use the local files.

We provide the script [`get_refs.py`](get_refs.py) to download these reference files and optionally pull container images to the location of your choice.
If you have downloaded the full `scpca-nf` repository, this script is included in the base directory.
Alternatively, you can download and this script on its own to the location of your choice with the following commands:

<!-- TODO: Update to `main` before merging -->
```sh
wget https://raw.githubusercontent.com/AlexsLemonade/scpca-nf/development/get_refs.py
chmod +x get_refs.py
```


Once you have downloaded the script and made it executable with the `chmod` command, running the script will download the files required for mapping gene expression data sets to the subdirectory `scpca-references` at your current location.
The script will also create a parameter file named `localref_params.yaml` that defines the `ref_rootdir` and `assembly` Nextflow parameter variables required to use these local data files.
To run with these settings

```sh
./get_refs.py
```



You can then direct Nextflow to use the parameters stored in `localref_params.yaml` by using the `-params-file` argument in a command such as the following:

```sh
nextflow run AlexsLemonade/scpca-nf \
  -params-file localref_params.yaml \
  -config user_template.config \
  -profile cluster
```

Note that other configuration settings such as [profiles](#setting-up-a-profile-in-the-configuration-file), must still be set in the configuration file directly.
However, you should **not** put `params.ref_rootdir` in the configuration file, as Nextflow may not properly create the sub-paths for the various reference files due to [Nextflow's precedence rules of setting parameters](https://www.nextflow.io/docs/latest/config.html#configuration-file).
The `ref_rootdir` parameter should *only* be specified in a parameter file or at the command line with the `--ref_rootdir` argument.

#### Additional reference files

If you wil be performing genetic demultiplexing for hashed samples, you will need STAR index files as well as the ones included by default.
To obtain these files, you can add the `--star_index` flag:

```sh
./get_refs.py --star_index
```

If you will be analyzing spatial expression data, you will also need the Cell Ranger index as well, which can be obtained by adding the `--cellranger_index` flag.

#### Downloading container images

If your compute nodes do not have internet access, you will likely have to pre-pull the required container images as well.
When doing this, it is important to be sure that you also specify the revision (version tag) of the `scpca-nf` workflow that you are using.
For example, if you would run `nextflow run AlexsLemonade/scpca-nf -r v0.4.0`, then you will want to set `-r v0.4.0` for `get_refs.py` as well to be sure you have the correct containers.
Be default,  `get_refs.py` will download files and images associated with the latest release.

If your system uses Docker, you can add the `--docker` flag:

```sh
./get_refs.py --docker
```

For Singularity, you can similarly use the `--singularity` flag to pull images and cache them for use by Nextflow.
These images will be placed by default in a `singularity` directory at your current location.
If you would like to store them in a different location, use the `--singularity_dir` argument to specify that path.
The example below stores the image files in `$HOME/singularity`.
You will also need to set the `singularity.cacheDir` variable to match this location in your [configuration file profile](#setting-up-a-profile-in-the-configuration-file).

```sh
./get_refs.py --singularity --singularity_dir "$HOME/singularity"
```

## Repeating mapping steps

By default, `scpca-nf` is set up to skip the `salmon` mapping steps for any libraries in which the output files from the mapping step exist in the `checkpoints` folder of the output directory (i.e. the `.rad` files from `salmon alevin` and `quant.sf` files from `salmon quant`).
If the `salmon` version and transcriptome index are unchanged, this will save substantial processing time and cost, and avoids some of the sensitivity of the caching system used by `nextflow -resume`, which can sometimes result in rerunning steps unnecessarily.
However, if there have been updates to the `scpca-nf` workflow that include changes to the salmon version or transcriptome index (or if you change those on your own), you may want to repeat the mapping process.

To force repeating the mapping process, use the `--repeat_mapping` flag at the command line:

```sh
nextflow run AlexsLemonade/scpca-nf \
  --repeat_mapping
```

## Special considerations for specific data types

### Libraries with additional feature data (CITE-seq or cellhash)

Libraries processed using multiple modalities, such as those that include runs of CITE-seq or cellhash tags, will require a file containing the barcode IDs and sequences.
The file location should be specified in the `feature_barcode_file` for each library as listed in the [metadata file](#prepare-the-metadata-file); multiple libraries can and should use the same `feature_barcode_file` if the same feature barcode sequences are expected.

The `feature_barcode_file` itself is a tab separated file with one line per barcode and no header.
The first column will contain the barcode or antibody ID and the second column the barcode nucleotide sequence.
For example:

```
TAG01	CATGTGAGCT
TAG02	TGTGAGGGTG
```


### Multiplexed (cellhash) libraries

When processing multiplexed libraries that combine multiple samples into a pooled single-cell or single-nuclei library, we perform both cellhash-based demultiplexing and genetic demultiplexing.

To support both of these demultiplexing strategies, we currently require *ALL* of the following for multiplexed libraries:

- A single-cell RNA-seq run of the pooled samples
- A matched cellhash sequencing run for the pooled samples
- Separate bulk RNA-seq libraries for each sample in the pool
- A TSV file, `feature_barcode_file`, defining the cellhash barcode sequences.
- A TSV file, `cellhash_pool_file` that defines the sample-barcode relationship for each library/pool of samples

The `feature_barcode_file` for each library should be listed in the [metadata file](#prepare-the-metadata-file).

The  `cellhash_pool_file` location will be defined as a parameter in the [configuration file](#configuration-files), and should contain information for all libraries to be processed.
This file will contain one row for each library-sample pair (i.e. a library containing 4 samples will have 4 rows, one for each sample within), and should contain the following required columns:

| column_id       | contents                                                       |
|-----------------|----------------------------------------------------------------|
| `scpca_library_id`| Multiplexed library ID matching values in the metadata file. |
| `scpca_sample_id` | Sample ID for a sample contained in the listed multiplexed library |
| `barcode_id`      | The barcode ID used for the sample within the library, as defined in `feature_barcode_file` |

Other columns may be included for reference (such as the `feature_barcode_file` associated with the library), but these will not be used directly.

We have provided an example multiplex pool file for reference that can be found in [`examples/example_multiplex_pools.tsv`](examples/example_multiplex_pools.tsv).

### Spatial transcriptomics libraries

To process spatial transcriptomic libraries, all FASTQ files for each sequencing run and the associated `.jpg` file must be inside the `files_directory` listed in the [metadata file](#prepare-the-metadata-file).
The metadata file must also contain columns with the `slide_section` and `slide_serial_number`.

You will also need to provide a [docker image](https://docs.docker.com/get-started/) that contains the [Space Ranger software from 10X Genomics](https://support.10xgenomics.com/spatial-gene-expression/software/downloads/latest).
For licensing reasons, we cannot provide a Docker container with Space Ranger for you.
As an example, the Dockerfile that we used to build Space Ranger can be found [here](https://github.com/AlexsLemonade/alsf-scpca/tree/main/images/spaceranger).

After building the docker image, you will need to push it to a [private docker registry](https://www.docker.com/blog/how-to-use-your-own-registry/) and set `params.SPACERANGER_CONTAINER` to the registry location and image id in the `user_template.config` file.
*Note: The workflow is currently set up to work only with spatial transcriptomic libraries produced from the [Visium Spatial Gene Expression protocol](https://www.10xgenomics.com/products/spatial-gene-expression) and has not been tested using output from other spatial transcriptomics methods.*



## Output files

Upon completion of the `scpca-nf` workflow, the results will be published to the specified `outdir`.
Within the `outdir`, two folders will be present, `results` and `checkpoints`.

The `results` folder will contain the final output files produced by the workflow and the files that are typically available for download on the ScPCA portal.

Within the `results` folder, all files pertaining to a specific sample will be nested within a folder labeled with the sample ID.
All files in that folder will be prefixed by library ID, with the following suffixes:  `_unfiltered.rds`, `_filtered.rds`, `_processed.rds`, `_metadata.json`, and `_qc.html`.
The `_unfiltered.rds`, `_filtered.rds`, and `_processed.rds` files contain the quantified gene expression data as a [`SingleCellExperiment` object](https://bioconductor.org/packages/release/bioc/html/SingleCellExperiment.html).
For more information on the contents of these files, see the [ScPCA portal docs section on single cell gene expression file contents](https://scpca.readthedocs.io/en/latest/sce_file_contents.html).

See below for the expected structure of the `results` folder:

```
results
└── sample_id
    ├── library_id_filtered.rds
    ├── library_id_metadata.json
    ├── library_id_processed.rds
    ├── library_id_qc.html
    └── library_id_unfiltered.rds
```

If bulk libraries were processed, a `bulk_quant.tsv` and `bulk_metadata.tsv` summarizing the counts data and metadata across all libraries will also be present in the `results` directory.

The `checkpoints` folder will contain intermediate files that are produced by individual steps of the workflow, including mapping with `salmon`.
The contents of this folder are used to allow restarting the workflow from internal checkpoints (in particular so the initial read mapping does not need to be repeated, see [repeating mapping steps](#repeating-mapping-steps)), and may contain log files and other outputs useful for troubleshooting or alternative analysis.

The `rad` folder (nested inside the `checkpoints` folder) contains the output from running [`salmon alevin`](https://salmon.readthedocs.io/en/latest/alevin.html) with the `--rad` flag.
If bulk libraries are processed, there will be an additional `salmon` folder that contains the output from running [`salmon quant`](https://salmon.readthedocs.io/en/latest/file_formats.html) on each library processed.

All files pertaining to a specific library will be nested within a folder labeled with the library ID.
Additionally, for each run, all files related to that run will be inside a folder labeled with the run ID followed by the type of run (i.e. `rna` or `features` for CITE-seq) and nested within the library ID folder.

See below for the expected structure of the `checkpoints` folder:

```
checkpoints
├── rad
│   ├── library01
│   │   ├── run01-rna
│   │   └── run02-features
│   └── library02
│       └── run03-rna
└── salmon
```

By default, the direct output from running [`alevin-fry`](https://alevin-fry.readthedocs.io/en/latest/index.html) is not provided.
Within `scpca-nf`, the [counts matrix output from `alevin-fry`](https://alevin-fry.readthedocs.io/en/latest/quant.html#output) is directly imported into R as a `SingleCellExperiment` object and can be obtained in the `_unfiltered.RDS` file.
If you would like to obtain all files typically output from running `alevin-fry`, you may run the workflow with the `--publish_fry_outs` option at the command line.
This will tell the workflow to save the `alevin-fry` outputs to a folder labeled `alevinfry` nested inside the `checkpoints` folder.

```sh
nextflow run AlexsLemonade/scpca-nf \
  -config <path to config file>  \
  -profile <name of profile> \
  --publish_fry_outs
```

If genetic demultiplexing was performed, there will also be a folder called `vireo` with the output from running [vireo](https://vireosnp.readthedocs.io/en/latest/index.html) using genotypes identified from the bulk RNA-seq.
Note that we do not output the genotype calls themselves for each sample or cell, as these may contain identifying information.
