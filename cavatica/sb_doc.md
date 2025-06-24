# `scpca-nf`

The  `scpca-nf` workflow is used to process 10x single-cell data as part of the [Single-cell Pediatric Cancer Atlas (ScPCA) project](https://scpca.alexslemonade.org/).


The workflow processes fastq files from single-cell and single-nuclei RNA-seq samples using [alevin-fry](https://alevin-fry.readthedocs.io/en/latest/) to create gene by cell matrices.
The workflow outputs gene expression data in two formats: as [`SingleCellExperiment` objects](https://www.bioconductor.org/packages/release/bioc/html/SingleCellExperiment.html) and as [`AnnData` objects](https://anndata.readthedocs.io/en/latest/).
Reads from samples are aligned using selective alignment, to an index with transcripts corresponding to spliced cDNA and to intronic regions, denoted by alevin-fry as `splici`.
These matrices are filtered and additional processing is performed to calculate quality control statistics, create reduced-dimension transformations, assign cell types using both [`SingleR`](https://bioconductor.org/packages/release/bioc/html/SingleR.html) and [`CellAssign`](https://docs.scvi-tools.org/en/stable/user_guide/models/cellassign.html), and create output reports.
`scpca-nf` can also process libraries with ADT tags (e.g., CITE-seq), multiplexed libraries (e.g., cell hashing), bulk RNA-seq, and spatial transcriptomics samples.

Full documentation for the workflow is available in the [scpca-nf GitHub repository](https://github.com/AlexsLemonade/scpca-nf).

## Instructions for use on Cavatica

### Prerequisites

If you are reading this on Cavatica, you should have already created a Cavatica project to run the workflow, and copied the `scpca-nf` app to your project.
If you have not yet done this, create a project in Cavatica and within that project go to the Apps tab, click "Add apps > Projects" and copy the `scpca-nf` app from the `scpca-nf-test` project (you will need to have been granted access to this project).

Make sure that your project settings include "Allow network access", so that the workflow can download the required reference files and software dependencies. 
We also highly recommend that you enable spot instances.

Finally, on your local machine, you will need to [install `sbpack`](https://docs.sevenbridges.com/docs/add-nextflow-apps-through-the-command-line-interface-cli#about-sbpack) in order to be able to use the `sbmanifest` tool.
You will also need to have set up [credentials for the Seven Bridges CLI](https://docs.sevenbridges.com/docs/store-credentials-to-access-seven-bridges-client-applications-and-libraries).


### Required input files

- FASTQ files for single-cell or single-nuclei RNA-seq libraries, or other supported modalities.
  - These files should be uploaded to Cavatica in the project that you are planning to use for running the workflow.
  - Each run should be in a separate directory, named with a unique identifier corresponding to the `scpca_run_id` in the run metadata file.
  - See [File organization](https://github.com/AlexsLemonade/scpca-nf/blob/main/external-instructions.md#file-organization) for more details.
- A run metadata file (TSV)
  - The run metadata file contains information about each run and library required for running the workflow.
    - Required columns are: `scpca_run_id`, `scpca_library_id`, `scpca_sample_id`, `scpca_project_id`, `technology`, `assay_ontology_term_id`, `seq_unit`, `sample_reference`, and `files_directory`.
    - For a complete listing of the column requirements and contents see [Prepare the sample metadata file](https://github.com/AlexsLemonade/scpca-nf/blob/main/external-instructions.md#prepare-the-run-metadata-file)
    - Note that you will likely want to specify the `files_directory` as the relative path _within the project_ to the directory containing the FASTQ files for the run; this will be transformed to a Cavatica `vs://` url by the `sbmanifest` tool
- A sample metadata file (TSV)
  - This file contains information about the biological samples, such as diagnosis, age, sex, cell line, or other relevant information.
  - A `scpca_sample_id` column is required, which should match the `scpca_sample_id` values in the run metadata file.
  - The sample metadata file should be uploaded to your Cavatica project.
  - See [Prepare the sample metadata file](https://github.com/AlexsLemonade/scpca-nf/blob/main/external-instructions.md#prepare-the-sample-metadata-file) for more details.

#### Additional files for specific modalities

- For ADT (CITE-seq) libraries, you will also need to provide a feature barcode file
  - The feature barcode file is a tab separated file with one line per barcode and no header, and should be uploaded to your Cavatica project.
  - The path to the feature barcode file must be included in the run metadata file, in the `feature_barcode_file` column.
  - See [Libraries with additional feature data](https://github.com/AlexsLemonade/scpca-nf/blob/main/external-instructions.md#libraries-with-additional-feature-data-adt-or-cellhash) for more details.


### Adding Cavatica paths to metadata files

To run the workflow on Cavatica, you will need to transform any relative paths in the run metadata file to Cavatica `vs://` uris.
This can be done using the `sbmanifest` tool, following the [Seven Bridges instructions](https:/docs.sevenbridges.com/docs/using-sample-sheets-with-nextflow-apps)

Briefly, you will need to run a command like the one below, replacing `{user/project}` with your Cavatica username and project name, and `run_metadata.tsv` with the path to your run metadata file.
The command below will validate that all of the paths in the `files_directory` column exist in your Cavatica project, and will output a new file `run_metadata-cavatica.tsv` with the transformed paths.

```
sbmanifest \
  --projectid "{user/project}" \
  --sample-sheet run_metadata.tsv \
  --output run_metadata-cavatica.tsv \
  --columns files_directory \
  --validate
```

If your run metadata file includes paths within your project such as for the `feature_barcode_file`, add those column names to the `--columns` argument to transform those paths as well: e.g., `--columns files_directory feature_barcode_file`.

Once you have validated and transformed the paths in your run metadata file, upload the new file to your Cavatica project (optionally using `sbmanifest` and the additional `--upload` flag).

### Running the workflow

To run the workflow, you can use the Cavatica web interface to launch the `scpca-nf` app.

- Go to the Apps tab in your project, and click the button "Run" next to the `scpca-nf` app.
- In the "Inputs" section, you will need to select the "Run metadata file" and "Sample metadata file" that you prepared earlier.
- All other fields are optional, but the following are the ones most likely to be useful:
  - "Output directory" can be changed to a name of your choice. The default is `scpca_out`, and Cavatica will create a unique name if the named  directory already exists.
  - "Perform cell typing" is False by default, as cell typing is slow, but it can be enabled if desired.
  - "Project ID(s)" and "Run IDs" can be used to limit the runs to a subset of those in the run metadata file. Each takes a comma-separated list of values.

Once you have set your desired options, click the "Run" button at the top of the page to launch the workflow!
When it is complete, you should see a link to the output directory, which will include a `checkpoints` directory with some intermediate files and the `results` directory with the final output files.
Refer to the main `scpca-nf` documentation for [more information about the output files](https://github.com/AlexsLemonade/scpca-nf/blob/main/external-instructions.md#output-files)
