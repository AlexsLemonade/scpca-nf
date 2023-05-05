# References for `scpca-nf`

This folder contains files and scripts related to maintaining the reference files available for use with `scpca-nf`.

1. `ref-metadata.tsv`: Each row of this TSV file corresponds to a reference that is available for mapping with `scpca-nf`.
The columns included specify the `organism` (e.g., `Homo_sapiens`), `assembly`(e.g.,`GRCh38`), and `version`(e.g., `104`) of the `fasta` obtained from [Ensembl](https://www.ensembl.org/index.html) that was used to build the reference files.
This file is used as input to the `build-index.nf` workflow, which will create all required index files for `scpca-nf` for the listed organisms in the metadata file, provided the `fasta` and `gtf` files are stored in the proper location on S3.
See [instructions for adding additional organisms](#adding-additional-organisms) for more details.

2. `scpca-refs.json`: Each entry of this file contains a supported reference for mapping with `scpca-nf` and the name used to refer to that supported reference, e.g., `Homo_sapiens.GRCh38.104`.
For each supported reference, a list of all the reference files that are needed to run `scpca-nf` will be included.
This file is required as input to `scpca-nf`.

## Adding additional organisms

Follow the below steps to add support for additional references:

1. Download the desired `fasta` and `gtf` files for the organism of choice from `Ensembl`.
Add these to the `S3://scpca-references` bucket with the following directory structure:
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
