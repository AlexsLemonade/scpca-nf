<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->

**Table of Contents**

- [Running `scpca-nf` as a Data Lab staff member](#running-scpca-nf-as-a-data-lab-staff-member)
  - [Additional flags and parameters](#additional-flags-and-parameters)
  - [Testing the workflow](#testing-the-workflow)
  - [Running `scpca-nf` for ScPCA Portal release](#running-scpca-nf-for-scpca-portal-release)
  - [Processing example data](#processing-example-data)
- [Maintaining references for `scpca-nf`](#maintaining-references-for-scpca-nf)
  - [Adding additional organisms](#adding-additional-organisms)
  - [Adding additional cell type references](#adding-additional-cell-type-references)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

## Running `scpca-nf` as a Data Lab staff member

This section provides instructions for running the main workflow, found in [`main.nf`](main.nf).
Note that there are two other workflows: [`build-index.nf`](build-index.nf) for building reference indices (see #adding-additional-organisms), and [`build-celltype-ref.nf`](build-cellltype-ref.nf) for creating cell type annotation references.

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

### Additional flags and parameters

There are several flags and/or parameters which you may additionally wish to specify, as follows.

- Nextflow flags:
  - `-resume`: Resume workflow from most recent checkpoint
  - `-with-tower`: Use `Nextflow Tower` to monitor workflow (requires separate [Nextflow Tower registration](https://tower.nf/))
- Workflow parameters:
  - `--run_ids list,of,ids`: A custom comma-separated list of ids (run, library, or sample) for this run.
  - `--project list,of,project_ids`: A custom comma-separated list of project ids for this run
    [The default](config/profile_ccdl.config) run ids are `"SCPCR000001,SCPCS000101"`.
  - `--repeat_mapping`: Use this flag to repeat mapping, even if results already exist.
    - By default, the workflow checks whether each library has existing `alevin-fry` or `salmon` mapping results, and skips mapping for libraries with existing results.
      Using this flag will override that default behavior and repeat mapping even if the given library's results exist.
    - For more implementation details, please refer to the [external instructions](external-instructions.md#repeating-mapping-steps).
  - `--skip_genetic_demux`: Use this flag to skip genetic demultiplexing, which is turned on by default.
    - Genetic demultiplexing requires mapping of both bulk and single-cell data, followed by SNP calling and genetic demultiplexing, which can be quite time consuming.
    - When genetic demultiplexing is skipped, the workflow will still perform cellhash-based demultiplexing, if available for a given library.
  - `--repeat_genetic_demux`: Use this flag to repeat genetic demultiplexing, even if results already exist.
    - By default, the workflow checks whether each library has existing genetic demultiplexing results, and skips genetic demultiplexing for libraries with existing results.
      Using this flag will override that default behavior and repeat genetic demultiplexing even if the given library's results exist.
  - `--perform_celltyping`: Use this flag to perform cell type annotation, which is turned off by default.
  - `--repeat_celltyping`: Use this flag to repeat cell type annotation, even if results already exist.
    - By default, the workflow checks whether each library has existing cell type annotation results for `SingleR` and/or `CellAssign` (depending on references for that library).
      Using this flag will override that default behavior and repeat cell type annotation even if the given library's results exist.
    - This flag is _only considered_ if `--perform_celltyping` is also used.

Please refer to [`nextflow.config`](nextflow.config) and [other configuration files](config/) for other parameters which can be modified.

### Testing the workflow

To test the workflow with the `stub`, use the following command:

```
nextflow -log stub-run.log run main.nf -stub -profile stub
```

The run log will be saved to `stub-run.log`.
By default, the output directory (`outdir` parameter) is set as `test/output`.
You may wish to verify that output contents in `test/output/checkpoints` and `test/output/results` are as expected when running the `stub` workflow.

Please refer to our [`CONTRIBUTING.md`](CONTRIBUTING.md#stub-workflows) for more information on maintaining the stub workflow.

### Running `scpca-nf` for ScPCA Portal release

When running the workflow for a project or group of samples that is ready to be released on ScPCA portal, please use the tag for the latest release:

```
nextflow run AlexsLemonade/scpca-nf -r v0.6.1 -profile ccdl,batch --project SCPCP000000
```

### Processing example data

We provide an [example of the expected outputs](./examples/README.md#example-output) after running `scpca-nf` available for external users.
If there have been major updates to the directory structure or the contents of the output, the example data should be re-processed such that the example output we provide mimics the current expected output from `scpca-nf`.

First, please check the metadata files present in `s3://scpca-references/example-data` are up to date with changes in the workflow.
There should be both an `example_run_metadata.tsv` and `example_sample_metadata.tsv`.
The columns of these files should match the expected input columns of the workflow (see the section on preparing the [run metadata](./external-instructions.md#prepare-the-run-metadata-file) and [sample metadata](./external-instructions.md#prepare-the-sample-metadata-file)).

Additionally, the `example_run_metadata.tsv` should contain at least 1 row with `run01` in the `scpca_run_id` column and `s3://scpca-references/example-data/example_fastqs/run01` in the `files_directory` column.

Once you have confirmed that the metadata looks correct, the following commands should be used to run the workflow and process the example data:

```sh
nextflow pull AlexsLemonade/scpca-nf -r development

nextflow run AlexsLemonade/scpca-nf -r development \
  -profile ccdl,batch \
  --run_ids run01 \
  --run_metafile s3://scpca-references/example-data/example_run_metadata.tsv \
  --sample_metafile s3://scpca-references/example-data/example_sample_metadata.tsv \
  --outdir s3://scpca-references/example-data/scpca_out
```

After successful completion of the run, the `scpca_out` folder containing the outputs from `scpca-nf` should be zipped up and stored at the following location: `s3://scpca-references/example-data/scpca_out.zip`.
Make sure to adjust the settings to make the zip file publicly accessible.

## Maintaining references for `scpca-nf`

Inside the `references` folder are files and scripts related to maintaining the reference files available for use with `scpca-nf`.

1. `ref-metadata.tsv`: Each row of this TSV file corresponds to a reference that is available for mapping with `scpca-nf`.
   The columns included specify the `organism` (e.g., `Homo_sapiens`), `assembly`(e.g.,`GRCh38`), and `version`(e.g., `104`) of the `fasta` obtained from [Ensembl](https://www.ensembl.org/index.html) that was used to build the reference files.
   This file is used as input to the `build-index.nf` workflow, which will create all required index files for `scpca-nf` for the listed organisms in the metadata file, provided the `fasta` and `gtf` files are stored in the proper location on S3.
   See [instructions for adding additional organisms](#adding-additional-organisms) for more details.

2. `scpca-refs.json`: Each entry of this file contains a supported reference for mapping with `scpca-nf` and the name used to refer to that supported reference, e.g., `Homo_sapiens.GRCh38.104`.
   For each supported reference, a list of all the reference files that are needed to run `scpca-nf` will be included.
   This file is required as input to `scpca-nf`.

3. `celltype-reference-metadata.tsv`: Each row of this TSV file corresponds to a supported cell type reference available for cell type annotation using `add-celltypes.nf`.
   For all references, the following columns will be populated: `celltype_ref_name`, `celltype_ref_source` (e.g., `celldex`), supported `celltype_method` (e.g., `SingleR`).
   All references obtained from the `PanglaoDB` source also require an `organs` column containing the list of supported `PanglaoDB` organs to include when building the reference.
   This should be a comma-separated list of all organs to include.
   To find all possible organs, see the `organs` column of `PanglaoDB_markers_2020-03-27.tsv`.
   This file is required as input to the `build-celltype-ref.nf` workflow, which will create all required cell type references for performing cell type annotation from the main workflow.
   See [instructions for adding additional cell type references](#adding-additional-cell-type-references) for more details.

4. `PanglaoDB_markers_2020-03-27.tsv`: This file is used to build the cell type references from `PanglaoDB`.
   This file was obtained from clicking the `get tsv file` button on the [PanglaoDB Dataset page](https://panglaodb.se/markers.html?cell_type=%27choose%27) and replacing the date in the filename with a date in ISO8601 format.
   This file is required as input to the `build-celltype-ref.nf` workflow, which will create all required cell type references for the main workflow to use during cell type annotation.

### Adding additional organisms

Adding additional organisms is handled, in part, by the `build-index.nf` workflow.

Follow the below steps to add support for additional references:

1. Download the desired `fasta` and `gtf` files for the organism of choice from `Ensembl`.
   Add these to the `S3://scpca-references` bucket with the following directory structure, where the root directory here corresponds to the `organism` and the subdirectory corresponds to the `Ensembl` version:

```
homo_sapiens
└── ensembl-104
    ├── annotation
    │   └── Homo_sapiens.GRCh38.104.gtf.gz
    └── fasta
        ├── Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai
        └── Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
```

2. Add the `organism`, `assembly`, and `version` associated with the new reference to the `ref-metadata.tsv` file.
3. Generate an updated `scpca-refs.json` by running the script, `create-reference-json.R`, located in the `scripts` directory.
4. Generate the index files using `nextflow run build-index.nf -profile ccdl,batch` from the root directory of this repository.
5. Ensure that the new reference files are public and in the correct location on S3 (`s3://scpca-references`).

### Adding additional cell type references

Adding additional organisms is handled, in part, by the `build-celltype-ref.nf` workflow.

Follow the below steps to add support for additional cell type references.
We currently only support `celldex` and `PanglaoDB` for reference sources for `SingleR` and `CellAssign` cell type annotation, respectively.

1. Add the `celltype_ref_name`, `celltype_ref_source`, `celltype_method`, and `organs` (if applicable) for the new reference to `celltype-reference-metadata.tsv`.
2. Generate the new cell type references using `nextflow run build-celltype-ref.nf -profile ccdl,batch` from the root directory of this repository.
3. Ensure that the new reference files are public and in the correct location on S3.
   Reference file names will be automatically determined by the `build-celltype-ref.nf`.

- `SingleR` reference files, which are the full reference datasets from the `celldex` package, should be in `s3://scpca-references/celltype/singler_references` and named as `<singler_ref_name>_<source>_<version>.rds`.
  Corresponding "trained" model files for use in the cell type annotation workflow should be stored in `s3://scpca-references/celltype/singler_models`, named as `<singler_ref_name>_<source>_<version>_model.rds`.
  Note that the workflow parameter `singler_label_name` will determine which `celldex` dataset label is used for annotation; by default, we use `"label.ont"` (ontology labels).
  - `<singler_ref_name>` represents the reference dataset name.
  - `<source>` represents the reference dataset source. Currently only `celldex` is supported.
  - `<version>` represents the source (`celldex`) version, where periods are replaced with dashes (e.g. version `x.y.z` would be represented as `x-y-z`).
- `CellAssign` organ-specific reference gene matrices should be in `s3://scpca-references/celltype/cellassign_references` named as `<cellassign_ref_name>_<source>_<date>.tsv`.
  - `<cellassign_ref_name>` represents the reference dataset name.
    Names are established by the Data Lab as `<tissue/organ>-compartment` to represent a set of markers for a given tissue/organ.
    The specific tissues/organs used for each set of markers is described in [`celltype-reference-metadata.tsv`](references/celltype-reference-metadata.tsv).
  - `<source>` represents the reference dataset source. Currently only `PanglaoDB` is supported.
  - `<date>` represents the source (`PanglaoDB`) download date in ISO8601 format.
