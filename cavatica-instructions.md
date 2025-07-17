
# How to use `scpca-nf` on Cavatica

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**

- [Prerequisites](#prerequisites)
- [Required input files](#required-input-files)
    - [Additional files for specific options](#additional-files-for-specific-options)
  - [Adding Cavatica paths to metadata files](#adding-cavatica-paths-to-metadata-files)
- [Running the workflow](#running-the-workflow)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->


## Prerequisites

To run the `scpca-nf` workflow on Cavatica, you must first confirm that you have the ability to run Nextflow workflows on Cavatica.
Access to this feature is granted by Seven Bridges administrators, and you can request access by contacting the Seven Bridges support team.

Once you have access to Nextflow workflows, you will most likely want to create a new project for your analysis.
Make sure that your project settings include "Allow network access", so that the workflow can download the required reference files and software dependencies.
We also highly recommend that you enable spot instances.

Within that project go to the Apps tab, click "Add apps > Projects" and copy the `scpca-nf` app from the `scpca-nf-test` project.
At the moment, you will need to have been granted access to the `scpca-nf-test` project by an ALSF Data Lab member.

Finally, on your local machine, you will need to [install `sbpack`](https://docs.sevenbridges.com/docs/add-nextflow-apps-through-the-command-line-interface-cli#about-sbpack) in order to be able to use the `sbmanifest` tool.
You will also need to have set up [credentials for the Seven Bridges CLI](https://docs.sevenbridges.com/docs/store-credentials-to-access-seven-bridges-client-applications-and-libraries).


## Required input files

The input files required to run the workflow on Cavatica are similar to those required for running on your own infrastructure, with some additional steps to ensure that paths are correctly formatted for Cavatica.
For more detailed information about formatting, you should refer to [the main instructions](external-instructions.md).

- FASTQ files for single-cell or single-nuclei RNA-seq libraries, or other supported modalities.
  - These files should be uploaded to Cavatica in the project that you are planning to use for running the workflow.
  - Each run should be in a separate directory, named with a unique identifier corresponding to the `scpca_run_id` in the run metadata file.
  - See [File organization](external-instructions.md#file-organization) for more details.
- A run metadata file (TSV)
  - The run metadata file contains information about each run and library required for running the workflow.
    - Required columns are: `scpca_run_id`, `scpca_library_id`, `scpca_sample_id`, `scpca_project_id`, `technology`, `assay_ontology_term_id`, `seq_unit`, `sample_reference`, and `files_directory`.
    - For a complete listing of the column requirements and contents see [Prepare the sample metadata file](external-instructions.md#prepare-the-run-metadata-file)
    - Note that you will likely want to specify the `files_directory` as the relative path _within the Cavatica project_ to the directory containing the FASTQ files for the run; this will be transformed to a Cavatica `vs:///` url by the [`sbmanifest` tool](#adding-cavatica-paths-to-metadata-files).
- A sample metadata file (TSV)
  - This file contains information about the biological samples, such as diagnosis, age, sex, cell line, or other relevant information.
  - A `scpca_sample_id` column is required, which should match the `scpca_sample_id` values in the run metadata file.
  - The sample metadata file should be uploaded to your Cavatica project.
  - See [Prepare the sample metadata file](https://github.com/AlexsLemonade/scpca-nf/blob/main/external-instructions.md#prepare-the-sample-metadata-file) for more details.

#### Additional files for specific options

- For ADT (CITE-seq) libraries or samples multiplexed via cell hashing, you will also need to provide a feature barcode file.
  - The feature barcode file is a tab separated file with one line per barcode and no header, and should be uploaded to your Cavatica project.
  - The path to the feature barcode file must be included in the run metadata file, in the `feature_barcode_file` column.
  - For multiplexed libraries, you will also need a `cellhash_pool_file` that defines the sample-barcode relationship for each library/pool of samples.

  - See [Libraries with additional feature data](https://github.com/AlexsLemonade/scpca-nf/blob/main/external-instructions.md#libraries-with-additional-feature-data-adt-or-cellhash) for more details.
- If you plan to perform cell typing, you will also need to provide a project cell type metadata file.
  - See [Preparing the project cell type metadata file](external-instructions.md#preparing-the-project-cell-type-metadata-file) for more details.
- You may also have a file that defines cell types for individual cells in the library.
  - For more information on this see [Providing existing cell type labels in the main instructions](external-instructions.md#providing-existing-cell-type-labels).
  - The path to this file should be included in the run metadata file, in the `submitter_cell_types_file` column.


### Adding Cavatica paths to metadata files

To run the workflow on Cavatica, you will need to transform any relative paths in the run metadata file to Cavatica `vs://` uris.
This can be done using the `sbmanifest` tool, following the [Seven Bridges instructions](https:/docs.sevenbridges.com/docs/using-sample-sheets-with-nextflow-apps)

Briefly, you will need to run a command like the one below, replacing `{user/project}` with your Cavatica username and project name, and `run_metadata.tsv` with the path to your run metadata file.
The command below will output a new file `run_metadata-cavatica.tsv` with the transformed paths, which will now have the form `vs:///Projects/{project_uuid}/relative/path/`.

```shell
sbmanifest \
  --projectid "{user/project}" \
  --sample-sheet run_metadata.tsv \
  --output run_metadata-sbmanifest.tsv \
  --columns files_directory
```

If your run metadata file includes other paths within your project, such as for the `feature_barcode_file`, add those column names to the `--columns` argument to transform those paths as well: e.g., `--columns files_directory feature_barcode_file`.

Similarly, if you have a multiplexed library and a corresponding multiplex pools file with a `feature_barcode_file` column, you will need to run the `sbmanifest` command with that file listed for the `--sample-sheet` argument as well, to transform the paths in the `feature_barcode_file` column.

**⚠️ Caution:** If you have `NA` values for samples that do not have paths, the `sbmanifest` command will transform those to paths of the form `vs:///Projects/{project_uuid}/NA`, which will cause errors when running the workflow (and fails if you try to use the tool's internal `--validate` option).
To undo this modification, you can use a command like the following to replace those paths with `NA` in the output file:

```bash
sed 's|vs:///Projects/[^/]*/NA|NA|g' run_metadata-sbmanifest.tsv > run_metadata-cavatica.tsv
```

Following the addition of the appropriate `vs://` paths, upload the modified run metadata file and any other required input files to your Cavatica project.


## Running the workflow

To run the workflow, use the Cavatica web interface to launch the `scpca-nf` app.

- Go to the Apps tab in your project, and click the button "Run" next to the `scpca-nf` app.
- In the "Inputs" section, select the "Run metadata file" and "Sample metadata file" that you prepared earlier.
- All other fields are generally optional, but the following are the ones most likely to be useful:
  - "Output directory" can be changed to a name of your choice. The default is `scpca_out`, and Cavatica will create a unique name if the named directory already exists.
  - "Perform cell typing" is False by default, as cell typing is slow, but it can be enabled if desired.
    - If you want to perform cell typing, you will also need to select a "Project cell type metadata file". For more information on this file see [Preparing the project cell type metadata file](external-instructions.md#preparing-the-project-cell-type-metadata-file).
  - "Project ID(s)" and "Run IDs" can be used to limit the runs to a subset of those in the run metadata file. Each takes a comma-separated list of values.

Once you have set your desired options, click the "Run" button at the top of the page to launch the workflow!
When it is complete, you should see a link to the output directory, which will include a `checkpoints` directory with some intermediate files and the `results` directory with the final output files.
Refer to the main `scpca-nf` documentation for [more information about the output files](external-instructions.md#output-files)
