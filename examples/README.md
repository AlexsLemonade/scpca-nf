# Example Metadata & Configuration Files for `scpca-nf`

## Example files

This directory contains the following example files:

| File name | Description |
|-----------|-------------|
| `user_template.config` | An example [Nextflow configuration file](../external-instructions.md#configuration-files) for the `scpca-nf` workflow |
| `example_run_metadata.tsv` | An example [run metadata file](../external-instructions.md#prepare-the-run-metadata-file) for the `scpca-nf` workflow. |
| `example_sample_metadata.tsv`| An example [sample metadata file](../external-instructions.md#prepare-the-sample-metadata-file) for the `scpca-nf` workflow. |
| `example_multiplex_pools.tsv` | An example [multiplex pool file](../external-instructions.md#multiplexed-cellhash-libraries) for the `scpca-nf` workflow.|
| `example_project_celltype_metadata.tsv` | An example [project cell type metadata file](../external-instructions.md#preparing-the-cell-type-project-metadata-file) for performing optional cell type annotation in the `scpca-nf` workflow. |

These files provide examples of expected formatting and content, but note that the specific values in these files may not be applicable or sufficient for running `scpca-nf` directly on your system.

Note that `example_run_metadata.tsv` additionally refers to several example files, whose descriptions can be found as follows:

- The format for the `example_barcode_files/cite_barcodes.tsv` file is described in the [ADT processing section of `external-instructions.md`](../external-instructions.md#libraries-with-additional-feature-data-adt-or-cellhash).
- The format for the `example_barcode_files/cellhash_barcodes.tsv` file is described in the [multiplexed libraries section of `external-instructions.md`](../external-instructions.md#multiplexed-cellhash-libraries).
- The format for the `example_metadata_files/submitter_celltypes.tsv` file is described in the [cell type annotation section of `external-instructions.md`](../external-instructions.md#providing-existing-cell-type-labels).


## Testing your setup with example data

:warning: These instructions are only intended to be used to test accurate set up of a configuration file.
Before following these instructions, please ensure that you have already set up your own [configuration file](../external-instructions.md#configuration-files) and have [created and named a profile to use](../external-instructions.md#setting-up-a-profile-in-the-configuration-file).

You can test your configuration setup by performing a test run with the example data that we have provided.

We recommend using the example 10X dataset from a [human glioblastoma donor that was processed using the 10X Genomics' Next GEM Single Cell 3' Reagent Kits v3.1](https://www.10xgenomics.com/resources/datasets/2-k-sorted-cells-from-human-glioblastoma-multiforme-3-v-3-1-3-1-standard-6-0-0)(note: you may be prompted to provide an email and register upon navigating to the 10X downloads site).
The fastq files for this example data can be downloaded from the following link (**note:** These files will take approximately 10 GB of disk space upon download and expanding the tar file): [Brain_Tumor_3p_fastqs.tar](https://cf.10xgenomics.com/samples/cell-exp/6.0.0/Brain_Tumor_3p/Brain_Tumor_3p_fastqs.tar).

Following download and unzipping of the fastq files, you will need to create a tab-separated values **run** metadata file that looks like the following:

| scpca_run_id | scpca_library_id | scpca_sample_id | scpca_project_id | technology | assay_ontology_term_id | seq_unit | sample_reference        | files_directory              |
| ------------ | ---------------- | --------------- | ---------------- | ---------- | ---------------------- | -------- | ----------------------- | ---------------------------- |
| run01        | library01        | sample01        | project01        | 10Xv3.1    | EFO:XXX                | cell     | Homo_sapiens.GRCh38.104 | /path/to/example_fastq_files |

Be sure to enter the **full** path to the directory containing the fastq files in the `files_directory` column.

You will also need to create a tab-separated values **sample** metadata file.
At a minimum, the sample metadata file must contain a column with `scpca_sample_id` as the header.
The contents of this column should contain all unique sample ids that are present in the `scpca_sample_id` column of the run metadata file.
Additionally, you may include columns `is_cell_line` and `is_xenograft` to indicate the sample type:

- `is_cell_line`: Use `TRUE` if the sample is from a cell line and `FALSE` otherwise.
Cell type annotation will not be performed for samples that are `TRUE`.
- `is_xenograft`: Use `TRUE` if the sample is from a patient-derived xenograft and `FALSE` otherwise.

This information will be reflected in the summary QC report.

Below is an example of a sample metadata file:

| scpca_sample_id |is_cell_line | is_xenograft | diagnosis    | age |
| --------------- |-------------| -------------|------------ | --- |
| sample01        | FALSE       | FALSE        | glioblastoma | 71  |

**Note that the `diagnosis` and `age` columns are shown as example sample metadata one might include in the sample metadata file.
The metadata file that you create does not need to match this exactly, but it must contain the required `scpca_sample_id` column.**

The following command can then be used to test your configuration setup with the example data:

```
nextflow run AlexsLemonade/scpca-nf \
  -config <path to config file>  \
  -profile <name of profile> \
  --run_metafile <path to run metadata file> \
  --sample_metafile <path to sample metadata file>
```

Where `<path to config file>` is the **relative** path to the configuration file that you have setup after following the instructions on [creating a configuration file](../external-instructions.md#configuration-files), `<name of profile>` is the name of the profile that you chose when creating a profile, `<path to run metadata file>` is the **full** path to the run metadata TSV you created, and `<path to sample metadata file>` is the **full** path to the sample metadata TSV you created.
For the [example configuration file that we provided](./user_template.config), we used the profile name `cluster` and would indicate that we would like to use that profile at the command line with `-profile cluster`.
For more detailed information on setting up the metadata file for your own data, see instructions on [preparing the run metadata file](../external-instructions.md#prepare-the-run-metadata-file) and [preparing the sample metadata file](../external-instructions.md/#prepare-the-sample-metadata-file).

## Example output

You can download an example of the expected output files here: [`scpca_out.zip`](https://s3.amazonaws.com/scpca-references/example-data/scpca_out.zip).
For more information on the file structure and what to expect see the description of the [output files](../external-instructions.md#output-files).
