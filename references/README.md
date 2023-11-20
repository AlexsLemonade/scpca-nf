This directory contains information and metadata about references used in the `scpca-nf` workflow.

If you are a Data Lab member looking for more information about maintaining or adding additional references, please refer to [`internal-instructions.md`](../internal-instructions.md).

### Transcriptome reference files

These files are used by `build-index.nf` to generate transcriptome reference files.

- `ref-metadata.tsv` provides transcriptome reference metadata for organisms considered in the `scpca-nf` workflow.
This file is used by `build-infex.nf` to generate transcriptome references.
- `scpca-refs.json` is produced by `scripts/create-reference-json.R` and provides relative file paths (within `s3://scpca-references`) for where references used in the `scpca-nf` workflow should be stored.


### Cell type reference files

These files are used by `build-celltype-ref.nf` to generate cell type annotation reference files.

- `celltype-reference-metadata.tsv` provides information about the name, source, and contents of cell type annotation references established by the Data Lab for use with `SingleR` and `CellAssign`.
- `PanglaoDB_markers_<date>` is the full set of marker genes exported from `PanglaoDB`, versioned at the given date.