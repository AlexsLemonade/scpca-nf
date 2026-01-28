<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**

- [Running `scpca-nf` as a Data Lab staff member](#running-scpca-nf-as-a-data-lab-staff-member)
- [Testing the workflow during development](#testing-the-workflow-during-development)
  - [Stub processes](#stub-processes)
  - [Test data runs](#test-data-runs)
- [Processing example data](#processing-example-data)
  - [Processing example 10x Flex data](#processing-example-10x-flex-data)
- [Maintaining references for `scpca-nf`](#maintaining-references-for-scpca-nf)
  - [Adding and updating reference files](#adding-and-updating-reference-files)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

## Running `scpca-nf` as a Data Lab staff member

The Data Lab will usually launch the `scpca-nf` workflows via GitHub Actions in the ScPCA-admin repository.
This is used to ensure that appropriate output locations and permissions are used.
For details about those runs, see the [Data Processing README](https://github.com/AlexsLemonade/ScPCA-admin/blob/main/data-processing-and-handoff/README.md)

Note that you may want to take advantage of the "Additional Nextflow options` field to set any parameters that are not already integrated with the GitHub Actions there.


## Testing the workflow during development

### Stub processes

Stub testing is useful for initial tests to be sure that the general workflow logic is correct.
However, it does not test the actual processing of data, so will miss errors in command execution or data processing.
To test the workflow with the `stub`, use the following command:

```
nextflow run main.nf -stub -profile stub -log stub-run.log
```

The run log will be saved to `stub-run.log`.
By default, the output directory (`outdir` parameter) is set as `test/output`.
You may wish to verify that output contents in `test/output/checkpoints` and `test/output/results` are as expected when running the `stub` workflow.

Please refer to our [`CONTRIBUTING.md`](CONTRIBUTING.md#stub-workflows) for more information on maintaining the stub workflow.


### Test data runs

When you are ready to test the workflow with real data, you will want to use the `testing` mode in GitHub Actions.
This can be configured to run the workflow from any branch using a subset of ScPCA data.
For more information, see the [Data Processing README](https://github.com/AlexsLemonade/ScPCA-admin/blob/main/data-processing-and-handoff/README.md
)

## Processing example data

We provide an [example of the expected outputs](./examples/README.md#example-output) after running `scpca-nf` available for external users.
If there have been major updates to the directory structure or the contents of the output, the example data should be re-processed such that the example output we provide mimics the current expected output from `scpca-nf`.

First, please check the metadata files present in `s3://scpca-nf-references/example-data` are up to date with changes in the workflow and reflect the contents of the files present in the `examples` directory of this repository.
Each of these files should be present, with the expected input columns as described in each documentation link.


- `example_run_metadata.tsv` ([documentation](./external-instructions.md#prepare-the-run-metadata-file))
- `example_sample_metadata.tsv` ([documentation](./external-instructions.md#prepare-the-sample-metadata-file))
- `example_project_celltype_metadata.tsv` ([documentation](./external-instructions.md#preparing-the-project-cell-type-metadata-file))
- `example_multiplex_pools.tsv` ([documentation](./external-instructions.md#multiplexed-cellhash-libraries))

Once you have confirmed that the metadata looks correct, use the [Run scpca-nf workflow on AWS Batch](https://github.com/AlexsLemonade/ScPCA-admin/actions/workflows/run-scpca-nf.yaml) GHA with the run mode set as `example` to process the example data.

### Processing example 10x Flex data

The example data described above includes both 10x Chromium and 10x Flex data.

There are two example datasets available on S3 that can be used specifically for testing changes to the `cellranger-flex.nf` module.
FASTQ files were downloaded from 10x Genomics, unzipped, and then copied to `s3://scpca-nf-references/example-data/example_fastqs`.
The information for these datasets were then added to `examples/example_run_metadata.tsv`, `examples/example_sample_metadata.tsv`, and `example_multiplex_pools.tsv`.
The datasets used are listed below:

1. [library06 - Human Kidney Nuclei - Singleplexed](https://10x.vercel.app/datasets/Human_Kidney_4k_GEM-X_Flex)
2. [library07 - Human PBMCs - Multiplexed](https://10x.vercel.app/datasets/80k_Human_PBMCs_PTG_MultiproPanel_IC_4plex)

For the second dataset (Human PBMCs), only the GEX FASTQ files were saved to S3.

To process only these datasets, use the `example` run mode with the [Run scpca-nf workflow on AWS Batch](https://github.com/AlexsLemonade/ScPCA-admin/actions/workflows/run-scpca-nf.yaml) GHA and specify the appropriate run ids, currently "library06,library07".

## Maintaining references for `scpca-nf`

Inside the `references` folder are files and scripts related to maintaining the reference files available for use with `scpca-nf`.

1. `ref-metadata.tsv`: Each row of this TSV file corresponds to a reference that is available for mapping with `scpca-nf`.
   The columns included specify the `organism` (e.g., `Homo_sapiens`), `assembly`(e.g.,`GRCh38`), and `version`(e.g., `104`) of the `fasta` obtained from [Ensembl](https://www.ensembl.org/index.html) that was used to build the reference files.
   This file is used as input to the `build-index.nf` workflow, which will create all required index files for `scpca-nf` for the listed organisms in the metadata file, provided the `fasta`, `gtf`, and `cytoband` (used to build `inferCNV` gene order files) files are stored in the proper location on S3.

2. `scpca-refs.json`: Each entry of this file contains a supported reference for mapping with `scpca-nf` and the name used to refer to that supported reference, e.g., `Homo_sapiens.GRCh38.104`.
   For each supported reference, a list of all the reference files that are needed to run `scpca-nf` will be included.
   This file is required as input to `scpca-nf`.

3. `celltype-reference-metadata.tsv`: Each row of this TSV file corresponds to a supported cell type reference available for cell type annotation.
   This file is required as input to the `build-celltype-ref.nf` workflow to create and/or update cell type references.
   For all references, the following columns must be populated:

     - `celltype_ref_name` (e.g., `BlueprintEncodeData` or `blood-compartment`)
     - `celltype_ref_source` (e.g., `celldex` or `PanglaoDB`)
     - `celltype_method` (e.g., `SingleR` or `CellAssign`)
     - All references obtained from the `PanglaoDB` source also require an `organs` column containing the list of supported `PanglaoDB` organs to include when building the reference.
       This should be a comma-separated list of all organs to include.
       To find all possible organs, see the `organs` column of `PanglaoDB_markers_2020-03-27.tsv`.

4. `PanglaoDB_markers_2020-03-27.tsv`: This file is used to build the cell type references from `PanglaoDB`.
   This file was obtained from clicking the `get tsv file` button on the [PanglaoDB Dataset page](https://panglaodb.se/markers.html?cell_type=%27choose%27) and replacing the date in the filename with a date in ISO8601 format.
   This file is required as input to the `build-celltype-ref.nf` workflow, which will create all required cell type references for the main workflow to use during cell type annotation.

5. The following files were generated in the `OpenScPCA-analysis` repository and copied to this repository for use in the workflow.
   They were initially obtained from the `OpenScPCA-analysis` repository at tag `v0.2.3`.

   The following files are used for cell typing, including assigning consensus cell types.
   If new cell typing methods are added or there are changes to references used for cell typing, these files will need to be updated.

   - [`panglao-cell-type-ontologies.tsv`](https://github.com/AlexsLemonade/OpenScPCA-analysis/blob/v0.2.3/analyses/cell-type-consensus/references/panglao-cell-type-ontologies.tsv)
   - [`scimilarity-mapped-ontologies.tsv`](https://github.com/AlexsLemonade/OpenScPCA-analysis/blob/v0.2.3/analyses/cell-type-scimilarity/references/scimilarity-mapped-ontologies.tsv)
   - [`consensus-cell-type-reference.tsv`](https://github.com/AlexsLemonade/OpenScPCA-analysis/blob/v0.2.3/analyses/cell-type-consensus/references/consensus-cell-type-reference.tsv)
   - [`consensus-validation-groups.tsv`](https://github.com/AlexsLemonade/OpenScPCA-analysis/blob/v0.2.3/analyses/cell-type-consensus/references/consensus-validation-groups.tsv)
   - [`validation-markers.tsv`](https://github.com/AlexsLemonade/OpenScPCA-analysis/blob/v0.2.3/analyses/cell-type-consensus/references/validation-markers.tsv)

   The following are used for `inferCNV` inference to determine which cell types to include in the normal reference.
   Additional rows will need to be added to these files if additional diagnoses are added to ScPCA.

   - [`diagnosis-groups.tsv`](https://github.com/AlexsLemonade/OpenScPCA-analysis/blob/v0.2.3/analyses/infercnv-consensus-cell-type/references/broad-diagnosis-map.tsv)
     - Before including in this repository, this file was modified from its `OpenScPCA-analysis` version to contain only two columns, `diagnosis_group` and `sample_diagnosis` (renamed from the original column name `submitted_diagnosis`), which are the column names that `scpca-nf` expects in this file.
   - [`diagnosis-celltype-groups.tsv`](https://github.com/AlexsLemonade/OpenScPCA-analysis/blob/v0.2.3/analyses/infercnv-consensus-cell-type/references/diagnosis-celltype-groups.tsv)

### Adding and updating reference files

Updates to the reference files in the `scpca-nf-references` bucket should be accomplished via the appropriate GitHub Actions in the ScPCA-admin repository.
For details, see https://github.com/AlexsLemonade/ScPCA-admin/blob/main/data-processing-and-handoff/README.md#updating-reference-files
