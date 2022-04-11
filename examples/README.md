# Example Metadata & Configuration Files for `scpca-nf`

This directory contains an example [metadata file](../external-data-instructions.md#prepare-the-metadata-file) and [configuration file](../external-data-instructions.md#configuration-files) for the `scpca-nf` workflow. 
These files should be used as an example of formats and content, but note that the values in these files may not be applicable or sufficient to allow running `scpca-nf` to be used directly on your system.

## Testing your setup with example data

If you would like to test your configuration setup, after following the instructions outlined in the [external instructions](../external-data-instructions.md), you can do so by performing a test run with the example data that we have provided.

We have an example 10X dataset from a [human glioblastoma donor that was processed using the 10X Genomics' Next GEM Single Cell 3' Reagent Kits v3.1](https://www.10xgenomics.com/resources/datasets/2-k-sorted-cells-from-human-glioblastoma-multiforme-3-v-3-1-3-1-standard-6-0-0) that is available for download. 
The fastq files for this example data can be downloaded from the following link: [Brain_Tumor_3p_fastqs.tar](https://cf.10xgenomics.com/samples/cell-exp/6.0.0/Brain_Tumor_3p/Brain_Tumor_3p_fastqs.tar).

Following download and unzipping of the fastq files, you will need to create a metadata file that looks like the following: 

| scpca_run_id | scpca_library_id | scpca_sample_id | scpca_project_id | technology | seq_unit | files_directory | 
| ------------ | ---------------- | --------------- | ---------------- | ---------- | -------- | --------------- | 
| run01 | library01 | sample01 | project01 | 10Xv3.1 | cell | /path/to/example_fastq_files | 

Be sure to enter the **full** path to the directory containing the fastq files in the `files_directory` column.

You can provide the path to the metadata file at the command line (here we name the file `example_metadata.tsv` but you can use the filename of your choice).
The following command can then be used to test your configuration setup with the example data: 

```
nextflow run AlexsLemonade/scpca-nf \
  -config my_config.config  \
  -profile <name of profile> \
  --run_metafile <path/to/example_metadata.tsv>
```

**Note:** This setup assumes that you have already created a configuration file, `my_config.config` (or the filename of your choice) and have created and named a profile to use. 
See the instructions on [creating a configuration file](../external-data-instructions.md#configuration-files) and [setting up a profile](../external-data-instructions.md#setting-up-a-profile-in-the-configuration-file) for more information.

## Example output 

You can download an example of the expected output files here: [scpca_out.zip](https://s3.amazonaws.com/scpca-references/example-data/scpca_out.zip). 
For more information on the file structure and what to expect see the description of the [output files](../external-data-instructions.md#output-files). 
