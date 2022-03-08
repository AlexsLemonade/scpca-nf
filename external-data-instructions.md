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

Using `scpca-nf` requires a metadata file where each sequencing run to be processed is a row and columns contain associated information about that run. 

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

To run the workflow, you will need to create a tab separated file with the following required columns: 

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
