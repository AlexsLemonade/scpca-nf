<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**

- [scpca-nf Data Lab Instructions](#scpca-nf-data-lab-instructions)
  - [Running scpca-nf as a Data Lab staff member](#running-scpca-nf-as-a-data-lab-staff-member)
    - [Processing example data](#processing-example-data)
  - [Maintaining references for `scpca-nf`](#maintaining-references-for-scpca-nf)
    - [Adding additional organisms](#adding-additional-organisms)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

# scpca-nf Data Lab Instructions

## Running scpca-nf as a Data Lab staff member

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

When running the workflow for a project or group of samples that is ready to be released on ScPCA portal, please use the tag for the latest release:

```
nextflow run AlexsLemonade/scpca-nf -r v0.5.3 -profile ccdl,batch --project SCPCP000000
```

### Processing example data

We provide an [example of the expected outputs](./examples/README.md#example-output) after running `scpca-nf` available for external users.
If there have been major updates to the directory structure or the contents of the output, the example data should be re-processed such that the example output we provide mimics the current expected output from `scpca-nf`.

The following command should be used to run the workflow and process the example data:

```sh
nextflow run AlexsLemonade/scpca-nf \
  -profile ccdl,batch \
  --run_ids run01 \
  --run_metafile s3://scpca-references/example-data/example_metadata.tsv \
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

3. `celltype-reference-metadata.tsv`: Each row of this TSV file corresponds to a supported cell type reference available for cell type assignment using `add-celltypes.nf`.
For all references, the following columns will be populated: `celltype_ref_name`, `celltype_ref_source` (e.g., `celldex`), supported `celltype_method` (e.g., `SingleR`).
All references obtained from the `PanglaoDB` source also require an `organs` column containing the list of supported `PanglaoDB` organs to include when building the reference.
This should be a comma-separated list of all organs to include.
To find all possible organs, see the `organs` column of `PanglaoDB_markers_27_Mar_2020.tsv`.
This file is required as input to the `build-celltype-ref.nf` workflow, which will create all required cell type references for `add-celltypes.nf`.
See [instructions for adding additional cell type references](#adding-additional-cell-type-references) for more details.

4. `PanglaoDB_markers_27_Mar_2020.tsv`: This file is used to build the cell type references from `PanglaoDB`.
This file was obtained from clicking the `get tsv file` button on the [PanglaoDB Dataset page](https://panglaodb.se/markers.html?cell_type=%27choose%27).
This file is required as input to the `build-celltype-ref.nf` workflow, which will create all required cell type references for `add-celltypes.nf`.

### Adding additional organisms

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

## Adding additional cell type references

Follow the below steps to add support for additional cell type references.
We currently only support `celldex` and `PanglaoDB` for reference sources.

1. Add the `celltype_ref_name`, `celltype_ref_source`, `celltype_method`, and `organs` (if applicable) for the new reference to `celltype-reference-metadata.tsv`.
2. Generate the new cell type references using `nextflow run build-celltype-ref.nf -profile ccdl,batch` from the root directory of this repository.
3. Ensure that the new reference files are public and in the correct location on S3 (`s3://scpca-references/celltype`).
