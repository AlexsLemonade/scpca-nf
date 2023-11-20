This directory contains information and metadata about references used in the `scpca-nf` workflow.

If you are a Data Lab member looking for more information about maintaining or adding additional references, please refer to [`internal-instructions.md`](../internal-instructions.md).

### Transcriptome reference files

These files are used by `build-index.nf` to generate transcriptome reference files.

- `ref-metadata.tsv` provides transcriptome reference metadata for organisms considered in the `scpca-nf` workflow.
This file is used by `build-infex.nf` to generate transcriptome references.
- `scpca-refs.json` contains the relative paths of the default reference files for mapping and quantification provided by the Data Lab, and is used for the `scpca-nf` configuration parameter `ref_json`.
  The base location for the reference files is given by the configuration parameter `ref_rootdir`.
  This file is produced by `scripts/create-reference-json.R`.


### Cell type reference files

These files are used by `build-celltype-ref.nf` to generate cell type annotation reference files.

- `celltype-reference-metadata.tsv` provides information about the name, source, and contents of cell type annotation references established by the Data Lab for use with `SingleR` and `CellAssign`.
- `PanglaoDB_markers_<date>` is the full set of marker genes exported from `PanglaoDB`, versioned at the given date.
