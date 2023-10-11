# Example Metadata & Configuration Files for `scpca-nf`

## Example files

This directory contains an example [metadata file](../external-data-instructions.md#prepare-the-metadata-file) and [configuration file](../external-data-instructions.md#configuration-files) for the `scpca-nf` workflow.
These files should be used as an example of formats and content, but note that the values in these files may not be applicable or sufficient to allow running `scpca-nf` to be used directly on your system.

## Testing your setup with example data

:warning: These instructions are only intended to be used to test accurate set up of a configuration file.
Before following these instructions, please ensure that you have already set up your own [configuration file](../external-data-instructions.md#configuration-files) and have [created and named a profile to use](../external-data-instructions.md#setting-up-a-profile-in-the-configuration-file).

You can test your configuration setup by performing a test run with the example data that we have provided.

We recommend using the example 10X dataset from a [human glioblastoma donor that was processed using the 10X Genomics' Next GEM Single Cell 3' Reagent Kits v3.1](https://www.10xgenomics.com/resources/datasets/2-k-sorted-cells-from-human-glioblastoma-multiforme-3-v-3-1-3-1-standard-6-0-0)(note: you may be prompted to provide an email and register upon navigating to the 10X downloads site).
The fastq files for this example data can be downloaded from the following link (**note:** These files will take approximately 10 GB of disk space upon download and expanding the tar file): [Brain_Tumor_3p_fastqs.tar](https://cf.10xgenomics.com/samples/cell-exp/6.0.0/Brain_Tumor_3p/Brain_Tumor_3p_fastqs.tar).


Following download and unzipping of the fastq files, you will need to create a tab-separated values **run** metadata file that looks like the following:

| scpca_run_id | scpca_library_id | scpca_sample_id | scpca_project_id | technology | assay_ontology_term_id | seq_unit | sample_reference | files_directory | submitter_cell_types_file |
| ------------ | ---------------- | --------------- | ---------------- | ---------- | ---------------------- | -------- | ---------------- | --------------- | ------------------------ |
| run01 | library01 | sample01 | project01 | 10Xv3.1 | EFO:XXX | cell | Homo_sapiens.GRCh38.104 | /path/to/example_fastq_files | /path/to/annotated_cell_types_file

Be sure to enter the **full** path to the directory containing the fastq files in the `files_directory` column, and similarly the full path for the `submitter_cell_types_file` column.

You will also need to create a tab-separated values **sample** metadata file.
At a minimum, the sample metadata file must contain a column with `scpca_sample_id` as the header_.
The contents of this column should contain all unique sample ids that are present in the `scpca_sample_id` column of the run metadata file.

Below is an example of a sample metadata file:

| scpca_sample_id | diagnosis | age |
| --------------- | --------- | --- |
| sample01 | glioblastoma | 71 |

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

Where `<path to config file>` is the **relative** path to the configuration file that you have setup after following the instructions on [creating a configuration file](../external-data-instructions.md#configuration-files), `<name of profile>` is the name of the profile that you chose when creating a profile, `<path to run metadata file>` is the **full** path to the run metadata TSV you created, and `<path to sample metadata file>` is the **full** path to the sample metadata TSV you created.
For the [example configuration file that we provided](./user_template.config), we used the profile name `cluster` and would indicate that we would like to use that profile at the command line with `-profile cluster`.
For more detailed information on setting up the metadata file for your own data, see instructions on [preparing the run metadata file](../external-data-instructions.md#prepare-the-run-metadata-file) and [preparing the sample metadata file](../external-instructions.md/#prepare-the-sample-metadata-file).

## Example output

You can download an example of the expected output files here: [`scpca_out.zip`](https://s3.amazonaws.com/scpca-references/example-data/scpca_out.zip).
For more information on the file structure and what to expect see the description of the [output files](../external-data-instructions.md#output-files).
