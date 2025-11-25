# `scpca-nf`

The  `scpca-nf` workflow is used to process 10x single-cell data as part of the [Single-cell Pediatric Cancer Atlas (ScPCA) project](https://scpca.alexslemonade.org/).

The workflow processes fastq files from single-cell and single-nuclei RNA-seq samples using [alevin-fry](https://alevin-fry.readthedocs.io/en/latest/) to create gene by cell matrices.
The workflow outputs gene expression data in two formats: as [`SingleCellExperiment` objects](https://www.bioconductor.org/packages/release/bioc/html/SingleCellExperiment.html) and as [`AnnData` objects](https://anndata.readthedocs.io/en/latest/).
Reads from samples are aligned using selective alignment, to an index with transcripts corresponding to spliced cDNA and to intronic regions, denoted by alevin-fry as `splici`.
These matrices are filtered and additional processing is performed to calculate quality control statistics, create reduced-dimension transformations, assign cell types using both [`SingleR`](https://bioconductor.org/packages/release/bioc/html/SingleR.html) and [`CellAssign`](https://docs.scvi-tools.org/en/stable/user_guide/models/cellassign.html), and create output reports.
`scpca-nf` can also process libraries with ADT tags (e.g., CITE-seq), multiplexed libraries (e.g., cell hashing), bulk RNA-seq, and spatial transcriptomics samples.

Full documentation for the workflow is available in the [scpca-nf GitHub repository](https://github.com/AlexsLemonade/scpca-nf).
Please also refer to the additional [instructions for use on Cavatica](https://github.com/AlexsLemonade/scpca-nf/blob/main/cavatica-instructions.md).

## Steps to run the workflow on Cavatica

To run the `scpca-nf` workflow on Cavatica, you will need to follow the following steps:

1. Upload all required input files to your Cavatica project
2. Prepare a metadata files for runs and samples, as well as any additional files that may be required for your samples. See [Required input files](https://github.com/AlexsLemonade/scpca-nf/blob/main/cavatica-instructions.md#required-input-files) for more details.
3. Transform the paths in the run metadata file to Cavatica `vs://` paths using the `sbmanifest` tool (see [Adding Cavatica paths to metadata files](https://github.com/AlexsLemonade/scpca-nf/blob/main/cavatica-instructions.md#adding-cavatica-paths-to-metadata-files) for more details).
4. Set options and run the `scpca-nf` workflow on Cavatica using the web interface.
