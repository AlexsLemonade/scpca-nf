
## CNV inference

`scpca-nf` currently performs copy-number variation (CNV) inference using [`inferCNV`](https://github.com/broadinstitute/infercnv).

By default, CNV inference is not performed.
As part of its algorithm, `inferCNV` uses a set of normal reference cells to quantify CNV events across all provided cells.
`scpca-nf` uses results from its cell type annotation module to create this normal reference.
As such, to perform CNV inference, you must also [perform cell type annotation, as described above](#cell-type-annotation).
Note that if CNV inference is turned on, cell type annotation will automatically be turned on as well.

You can turn on CNV inference by taking the following steps:

1. Ensure your [prepared sample metadata file](#prepare-the-sample-metadata-file) contains a `diagnosis` column.
   This column is **required** in your sample metadata file as `scpca-nf` uses this information to create the normal reference.
2. [Determine which cell types](#determining-normal-reference-cell-types) should be used for the `inferCNV` normal reference.
3. Using the information from the previous steps, [prepare the diagnosis metadata files](#preparing-diagnosis-metadata-files) which the workflow will use to create the normal cell type reference.
   You will need to provide the path/uri to this files as a workflow parameters, which you will need to define in your configuration file.
   For more information on adding parameters to your configuration file, see [Configuring scpca-nf for your environment](#configuring-scpca-nf-for-your-environment).
4. Run the workflow with the `--perform_cnv_inference` flag.

Once you have followed both the above steps and the [steps to prepare for cell type annotation](#cell-type-annotation) and you have added the path/uri to the diagnosis metadata files to your configuration file, you can use the following command to run the workflow with CNV inference:

```sh
nextflow run AlexsLemonade/scpca-nf \
  --perform_cnv_inference
```

### Determining normal reference cell types

`scpca-nf` uses the consensus cell type labels derived from `SingleR` and `CellAssign` annotations to create the set of normal reference cells.
The specific consensus cell types to use are determined by the sample's diagnosis.

To manage many diagnoses and consensus cell types, the Data Lab has organized [individual sample diagnoses into broad diagnosis groups](references/broad-diagnosis-map.tsv), and similarly [consensus cell types into validation groups](references/consensus-validation-groups.tsv).

`scpca-nf` uses these grouping to specify the normal reference cell types as follows:

1. Given a sample's diagnosis, identify the corresponding broad diagnosis
    * This step utilizes the provided `params.diagnosis_groups_file` TSV
2. Identify the consensus cell type validation groups that are expected to be normal cells in samples of that diagnosis
    * This step utilizes the provided `params.diagnosis_celltypes_file` TSV
3. Identify the individual consensus cell types that comprise those validation groups

To support the first two steps of this procedure, you will need to provide two metadata TSV files:

* [The diagnosis groups metadata file](#diagnosis-groups-metadata-file), which maps individual sample diagnoses to broad diagnosis groups
* [The diagnosis cell types metadata file](#diagnosis-cell-types-metadata-file) that maps broad diagnosis groups to consensus cell type validation groups

Instructions for creating these files are below.

### Preparing the diagnosis metadata files

#### Diagnosis groups metadata file

This path/uri to this file should be specified using the workflow parameter `diagnosis_groups_file`.
This file contain these two columns with the following information:

| column_id             | contents |
| --------------------- | -------- |
| `diagnosis_group`    | Broad diagnosis group |
| `submitted_diagnosis` | A sample-level diagnosis provided in the [sample metadata file](#prepare-the-sample-metadata-file) |

Please ensure there is only one submitted diagnosis per row; you can add additional rows if/when multiple `submitted_diagnosis` values correspond to the same `diagnosis_group`

We have provided an example diagnosis group metadata file for reference.

| [View example `diagnosis_groups.tsv` file](examples/example_diagnosis_groups.tsv) |
| --------------------------------------------------------------------------------------------------- |


#### Diagnosis cell types metadata file

This path/uri to this file should be specified using the workflow parameter `diagnosis_celltypes_file`.
This file contain these two columns with the following information:

|  |   |
| --- | --- |
| `diagnosis_group` | Broad diagnosis group. This should have a corresponding value in the diagnosis group metadata file |
| `celltype_groups` | Consensus cell type validation groups to include in an `inferCNV` normal reference for samples of the given broad diagnosis group  |


In this file, there should be only one row per `diagnosis_group`.

The `celltype_groups` column should be provided as a comma-separated list without spaces.
You can find the possible values to include in this column from the reference file [`references/consensus-validation-groups.tsv`](references/consensus-validation-groups.tsv).
This file maps individual consensus cell types (column `consensus_annotation`) to validation groups (column `validation_group_annotation`).
Values from the `validation_group_annotation` column should be used to populate the `celltype_groups` column in your diagnosis celltypes metadata file.


We have provided an example diagnosis cell type metadata file for reference.

| [View example `diagnosis_celltypes.tsv` file](examples/example_diagnosis_celltypes.tsv) |
| --------------------------------------------------------------------------------------------------- |

### Repeating CNV inference

TODO: https://github.com/AlexsLemonade/scpca-nf/issues/1025
