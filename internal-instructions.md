<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**

- [Running `scpca-nf` as a Data Lab staff member](#running-scpca-nf-as-a-data-lab-staff-member)
  - [Additional flags and parameters](#additional-flags-and-parameters)
  - [Testing the workflow with stub processes](#testing-the-workflow-with-stub-processes)
  - [Running `scpca-nf` for ScPCA Portal release](#running-scpca-nf-for-scpca-portal-release)
  - [Processing example data](#processing-example-data)
    - [Processing example 10x Flex data](#processing-example-10x-flex-data)
- [Maintaining references for `scpca-nf`](#maintaining-references-for-scpca-nf)
  - [Adding additional organisms](#adding-additional-organisms)
  - [Adding additional cell type references](#adding-additional-cell-type-references)
  - [Adding additional gene order files](#adding-additional-gene-order-files)
- [Running the merge workflow](#running-the-merge-workflow)
- [Comparing results between `scpca-nf` runs](#comparing-results-between-scpca-nf-runs)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

## Running `scpca-nf` as a Data Lab staff member

This section provides instructions for running the main workflow, found in [`main.nf`](main.nf).
Note that there are three other workflows:

- [`build-index.nf`](build-index.nf) for building reference indices (see #adding-additional-organisms)
- [`build-celltype-ref.nf`](build-celltype-ref.nf) for creating cell type annotation references (see #adding-additional-cell-type-references)
- [`merge.nf`](merge.nf) for merging processed objects produced by the main workflow (see #running-the-merge-workflow)

The instructions below assume that you are a member of the Data Lab with access to AWS.
Most of the workflow settings described are configured for the ALSF Childhood Cancer Data Lab computational infrastructure.
To process samples that are not part of the ScPCA project, please see the [instructions on using `scpca-nf` with external data](external-instructions.md).

To process single-cell and single-nuclei samples using `scpca-nf` you will need access to at least 24 GB of RAM and 12 CPUs, so we recommend using AWS batch.
The first step in running the workflow is ensuring that your AWS credentials are configured.

You can then run the workflow with the `batch` profile, which has been named in the `nextflow.config` file with full settings set up in the `profile_awsbatch.config` file.
Note that you will also need to use the `ccdl` profile to set input and output paths, and you can specify both profiles by separating them with a comma.

You can then run the workflow use the following command:

```
nextflow run AlexsLemonade/scpca-nf -profile ccdl,batch
```

Note that not all projects and samples are available at the paths specified in the `ccdl` profile.
This profile is primarily intended for development and testing purposes.
See the [Running `scpca-nf` for ScPCA Portal release](#running-scpca-nf-for-scpca-portal-release) section for details on full runs.

### Additional flags and parameters

There are several flags and/or parameters which you may additionally wish to specify, as follows.

- Nextflow flags:
  - `-resume`: Resume workflow from most recent cached results
  - `-with-tower`: Use `Nextflow Tower` to monitor workflow (requires separate [Nextflow Tower registration](https://tower.nf/))
- Workflow parameters:
  - `--run_ids list,of,ids`: A custom comma-separated list of ids (run, library, or sample) for this run.
  - `--project list,of,project_ids`: A custom comma-separated list of project ids for this run
    [The default](config/ccdl_profiles.config) run ids are `"SCPCR000001,SCPCS000101"`.
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
  - `--perform_cnv_inference`: Use this flag to perform CNV inference, which is turned off by default.
    - If CNV inference is specified, cell type annotation will automatically be turned on (e.g., `perform_celltyping` will be set to `true`)
  - `--repeat_cnv_inference`: Use this flag to repeat CNV inference, even if results already exist.
    - By default, the workflow checks whether each library has existing `inferCNV` results.
      Using this flag will override that default behavior and repeat CNV inference even if the given library's results exist.
    - This flag is _only considered_ if `--perform_cnv_inference` is also used.

Please refer to [`nextflow.config`](nextflow.config) and [other configuration files](config/) for other parameters which can be modified.

### Testing the workflow with stub processes

Stub testing is useful for initial tests to be sure that the general workflow logic is correct.
However, it does not test the actual processing of data, so will miss errors in command execution or data processing.
To test the workflow with the `stub`, use the following command:

```
nextflow -log stub-run.log run main.nf -stub -profile stub
```

The run log will be saved to `stub-run.log`.
By default, the output directory (`outdir` parameter) is set as `test/output`.
You may wish to verify that output contents in `test/output/checkpoints` and `test/output/results` are as expected when running the `stub` workflow.

Please refer to our [`CONTRIBUTING.md`](CONTRIBUTING.md#stub-workflows) for more information on maintaining the stub workflow.

### Running `scpca-nf` for ScPCA Portal release

When preparing for a new portal release, the workflow should first be run using the `ccdl_staging` profile.
Be sure to use the `-r` flag to specify the latest release tag for the workflow, and use the `--project` flag to specify the project id(s) to process.

For example:

```
nextflow run AlexsLemonade/scpca-nf -r v0.8.8 -profile ccdl_staging,batch --project SCPCP000000
```

When that run has completed successfully, check that the outputs are as expected.
You can then run the workflow using the `ccdl_prod` profile:

```
nextflow run AlexsLemonade/scpca-nf -r v0.8.8 -profile ccdl_prod,batch --project SCPCP000000
```

Both of these profiles have `-with-tower` set by default, and will use the [ScPCA workspace](https://cloud.seqera.io/orgs/CCDL/workspaces/ScPCA/watch) for monitoring (allowing all team members to see progress).
*Be sure to have set the `TOWER_ACCESS_TOKEN` environment variable to your Seqera platform token when triggering a run, or these runs will fail.*
The best way to do this is to store the token in 1Password, set the `TOWER_ACCESS_TOKEN` environment variable to the appropriate [secret reference](https://developer.1password.com/docs/cli/secret-references/), and then run the Nextflow command as above but prefixed with `op run -- `.

### Processing example data

We provide an [example of the expected outputs](./examples/README.md#example-output) after running `scpca-nf` available for external users.
If there have been major updates to the directory structure or the contents of the output, the example data should be re-processed such that the example output we provide mimics the current expected output from `scpca-nf`.

First, please check the metadata files present in `s3://scpca-references/example-data` are up to date with changes in the workflow and reflect the contents of the files present in the `examples` directory of this repository.
Each of these files should be present, with the expected input columns as described in each documentation link.

- `example_run_metadata.tsv` ([documentation](./external-instructions.md#prepare-the-run-metadata-file))
- `example_sample_metadata.tsv` ([documentation](./external-instructions.md#prepare-the-sample-metadata-file))
- `example_project_celltype_metadata.tsv` ([documentation](./external-instructions.md#preparing-the-project-cell-type-metadata-file))
- `example_multiplex_pools.tsv` ([documentation](./external-instructions.md#multiplexed-cellhash-libraries))

Once you have confirmed that the metadata looks correct, use the following commands to run the workflow and re-process the example data:

```sh
# Obtain the latest development version
nextflow pull AlexsLemonade/scpca-nf -r development

# Run the workflow with the example config
nextflow run AlexsLemonade/scpca-nf -r development -profile example,batch
```

After successful completion of the run, the `scpca_out` folder containing the outputs from `scpca-nf` should be zipped up and stored at the following location: `s3://scpca-references/example-data/scpca_out.zip`.
Be sure that only the results from `run01`, which is from `library01` and `sample01`, are included in the folder.
Make sure to adjust the settings to make the zip file publicly accessible.

#### Processing example 10x Flex data

Any samples that are processed using the [GEM-X Flex Gene Expression protocol from 10x Genomics](https://www.10xgenomics.com/products/flex-gene-expression) are quantified using [`cellranger multi`](https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/running-pipelines/cr-flex-multi-frp) instead of `alevin-fry`.

There are two example datasets available on S3 that can be used for testing changes to the `cellranger-flex.nf` module.
FASTQ files were downloaded from 10x Genomics, unzipped, and then copied to `s3://scpca-references/example-data/example_fastqs`.
The information for these datasets were then added to `examples/example_run_metadata.tsv`, `examples/example_sample_metadata.tsv`, and `example_multiplex_pools.tsv`.
The datasets used are listed below:

1. [library06 - Human Kidney Nuclei - Singleplexed](https://10x.vercel.app/datasets/Human_Kidney_4k_GEM-X_Flex)
2. [library07 - Human PBMCs - Multiplexed](https://10x.vercel.app/datasets/80k_Human_PBMCs_PTG_MultiproPanel_IC_4plex)

For the second dataset (Human PBMCs), only the GEX FASTQ files were saved to S3.

To process these datasets, use the `example` profile and specify the appropriate run, library, or sample IDs:

```sh
nextflow run AlexsLemonade/scpca-nf -r {branch or revision} -profile example,batch --run_ids library06,library07
```

## Maintaining references for `scpca-nf`

Inside the `references` folder are files and scripts related to maintaining the reference files available for use with `scpca-nf`.

1. `ref-metadata.tsv`: Each row of this TSV file corresponds to a reference that is available for mapping with `scpca-nf`.
   The columns included specify the `organism` (e.g., `Homo_sapiens`), `assembly`(e.g.,`GRCh38`), and `version`(e.g., `104`) of the `fasta` obtained from [Ensembl](https://www.ensembl.org/index.html) that was used to build the reference files.
   This file is used as input to the `build-index.nf` workflow, which will create all required index files for `scpca-nf` for the listed organisms in the metadata file, provided the `fasta`, `gtf`, and `cytoband` (used to build `inferCNV` gene order files) files are stored in the proper location on S3.
   See [instructions for adding additional organisms](#adding-additional-organisms) for more details.

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

   See [instructions for adding additional cell type references](#adding-additional-cell-type-references) for more details.

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

### Adding additional organisms

Adding additional organisms is handled, in part, by the `build-index.nf` workflow.

Follow the below steps to add support for additional references:

1. Download the desired `fasta` and `gtf` files for the organism of choice from `Ensembl`.
   Add these to the `s3://scpca-references` bucket with the following directory structure, where the root directory here corresponds to the `organism` and the subdirectory corresponds to the `Ensembl` version:

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
Specify which indexes should be built for this reference version, using the `include_salmon`, `include_cellranger`, and `include_star` columns.
If the `inferCNV` gene order file is also needed, also set the `include_infercnv` column to `TRUE` and follow instructions in the [Adding additional gene order files section](#adding-additional-gene-order-files) below.
3. Generate an updated `scpca-refs.json` by running the script, `create-reference-json.R`, located in the `scripts` directory.
4. Generate the index files using `nextflow run build-index.nf -profile ccdl,batch` from the root directory of this repository.
To generate the index files for only the new organism, use the `--build_refs` argument at the command line and specify the name of the reference to build, e.g., `nextflow run build-index.nf -profile ccdl,batch --build_refs Homo_sapiens.GRCh38.104`.
5. Ensure that the new reference files are public and in the correct location on S3 (`s3://scpca-references`).

### Adding additional cell type references

Adding additional references to use for cell type annotation is handled by the `build-celltype-ref.nf` workflow.

Reference files are created and automatically named by the `build-celltype-ref.nf`.
We currently only support `celldex` and `PanglaoDB` for reference sources for `SingleR` and `CellAssign` cell type annotation, respectively.

Follow these steps to add support for additional cell type references.

1. Add the `celltype_ref_name`, `celltype_ref_source`, `celltype_method`, and `organs` (if applicable) for the new reference to [`celltype-reference-metadata.tsv`](references/celltype-reference-metadata.tsv).

    - `{celltype_ref_name}` represents the reference dataset name.
      For use with `SingleR`, this should be taken directly from a `celldex` dataset.
      For `CellAssign`, names are established by the Data Lab as `{tissue/organ}-compartment` to represent a set of markers for a given tissue/organ.
    - `{celltype_ref_source}` represents the reference dataset source. Currently only `celldex` and `PanglaoDB` are supported for `SingleR` and `CellAssign`, respectively.
    - `{celltype_method}` represents which annotation method to use with the specified reference, either `SingleR` or `CellAssign`.
    - `organs` indicates which organs to be included in creation of references with `PanglaoDB` as the `celltype_ref_source`.
       This must be a comma separated list of all organs to include.

2. Generate the new cell type reference using `nextflow run build-celltype-ref.nf -profile ccdl,batch` from the root directory of this repository.
3. Ensure that the new reference files are public and in the correct location on S3.

`SingleR` reference files, which are the full reference datasets from the `celldex` package, should be in `s3://scpca-references/celltype/singler_references` and named as `{celltype_ref_name}_{celltype_ref_source}_{version}.rds`.
Corresponding "trained" model files for use in the cell type annotation workflow should be stored in `s3://scpca-references/celltype/singler_models`, named as `{celltype_ref_name}_{celltype_ref_source}_{version}_{gene_set_version}_{date_generated}_model.rds`.

  - `{celltype_ref_name}` is a given `celldex` dataset.
    - Note that the workflow parameter `singler_label_name` will determine which `celldex` dataset label is used for annotation; by default, we use `"label.ont"` (ontology labels).
  - `{celltype_ref_source}` is `celldex`.
  - `{version}` is the `celldex` version used during reference building, where we use dashes in place of periods (e.g., version `x.y.z` would be represented as `x-y-z`).
  - `{gene_set_version}` refers to the reference transcriptome or probe set used for mapping.
Currently, one model for the transcriptome references and one model for the flex probe sets are saved.

`CellAssign` organ-specific reference gene matrices should be stored in `s3://scpca-references/celltype/cellassign_references` and named as `{celltype_ref_name}_{celltype_ref_source}_{date}.tsv`.

  - `{celltype_ref_name}` is a given reference name established by the Data Lab.
  - `{celltype_ref_source}` is `PanglaoDB`
  - `{date}` is the `PanglaoDB` date, which serves as their version, in ISO8601 format.

### Adding additional gene order files

Adding additional [gene order files](https://github.com/broadinstitute/inferCNV/wiki/File-Definitions#gene-ordering-file) to use with `inferCNV` is handled by the `build-ref.nf` workflow.
A new gene order file will be needed if new a `gtf` file (e.g. a new `Ensembl` version or species) is added to references.
The gene order files consider chromosome arms as well, so an appropriate `cytoband` file with chromosome arm boundaries is also needed.

Follow these steps to add additional gene order files.

1. Ensure input files needed for creating the new gene order file have been added to `s3://scpca-references`.
This includes both a `gtf` file which can be downloaded from `Ensembl`, and a `cytoband` file delimiting chromosome arms which can be downloaded from the appropriate reference subdirectory in `ftp://hgdownload.cse.ucsc.edu/goldenPath/`.
When adding these files to the `s3://scpca-references` bucket, the directory structure should be as follows, where the root directory here corresponds to the `organism` and the subdirectory corresponds to the `Ensembl` version:

```
homo_sapiens
└── ensembl-104
    └── annotation
        ├── Homo_sapiens.GRCh38.104.gtf.gz
        └── Homo_sapiens.GRCh38.104_cytoband.txt.gz
```

If the `gtf` file is also new, be sure to also follow [these previous instructions](#adding-additional-organisms) for adding additional organisms.

2. Update the `include_infercnv` field in the `ref-metadata.tsv` file to `TRUE` for the new version being generated.
3. Generate an updated `scpca-refs.json` by running the script, `create-reference-json.R`, located in the `scripts` directory.
4. Generate the gene order file using `nextflow run build-index.nf -profile ccdl,batch` from the root directory of this repository.
To generate the index files for only the new organism, use the `--build_refs` argument at the command line and specify the name of the reference to build, e.g., `nextflow run build-index.nf -profile ccdl,batch --build_refs Homo_sapiens.GRCh38.104`.
5. Ensure that the gene order file is public and in the correct location on S3 (`s3://scpca-references`), for example:
```
homo_sapiens
└── ensembl-104
    └── infercnv
        └── Homo_sapiens.GRCh38.104_gene_order_arms.txt.gz
```

## Running the merge workflow

The workflow in `merge.nf` merges (but **does not integrate**) processed objects produced by the main `scpca-nf` workflow.
For a given SCPCA project id, the workflow creates a merged `SCE` object, a merged `AnnData` object, and an associated merged object HTML report.
Note that the workflow will also merge any present alternative experiments (e.g., ADT data from CITE-seq) _except_ for cell hashing alternative experiments.

The merge workflow requires two parameters:

- `project`, the SCPCA project id whose objects should be merged.
If running multiple projects, provide a comma separated list of project IDs.
- `merge_run_ids`, the set of run ids, library ids, or sample ids to include in the merged object.
By default all runs ids for the given project will be included
- `run_metafile`, the metadata file (`scpca-library-metadata.tsv`) which contains information about libraries to merge
  - This is specified in the `ccdl` profile configuration file

Data Lab members with access to AWS can run the workflow with the following command(s):

```
# Run a single project
nextflow run merge.nf -profile ccdl,batch --project SCPCP000000

# Run more than one project
nextflow run merge.nf -profile ccdl,batch --project SCPCP00000X,SCPCP00000Y

# Specify a set of run ids to use
nextflow run merge.nf -profile ccdl,batch --project SCPCP000000 --run_ids SCPCR00000X,SCPCR00000Y
```

## Comparing results between `scpca-nf` runs

To facilitate comparisons between the production version of `scpca-nf` results and the results from a staging (or other) run, we have written a script and notebook to compare metrics between different runs, producing an HTML report.
This script is found in the `scripts/compare-metrics` directory.

The script can be run for a single project using a command like the following, assuming that AWS credentials are available in the environment:

```
Rscript scripts/compare-metrics/compare-metrics.R \
  --project_id SCPCP000001 \
  --output_file "SCPCP000001_metrics_comparison.html"
```


Multiple projects can be included using a comma-separated list for the `--project_id` argument, and if the `--project_id` argument is omitted, the comparison will include all projects.

The default comparison will be between the metrics files available in the production and staging directories, treating the production versions as the reference run.
The files are expected to be found at `s3://nextflow-ccdl-results/scpca-prod/results` and `s3://nextflow-ccdl-results/scpca-staging/results`, by default, but these can be changed using the `--ref_s3` and `--comp_s3` arguments.
