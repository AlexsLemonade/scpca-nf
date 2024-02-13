# Overview

The `scpca-nf` workflow includes the option to perform cell type annotation, [as described in these instructions](./external-instructions.md#cell-type-annotation), with the reference-based method [`SingleR`](https://bioconductor.org/packages/release/bioc/html/SingleR.html) and/or the marker-gene based method [`CellAssign`](https://github.com/Irrationone/cellassign).

Several [reference datasets](./external-instructions.md#singler-references) and [marker gene lists](./external-instructions.md#cellassign-references) have been pre-compiled for use.

This document provides instructions for creating and using your own reference data set for `SingleR` annotation or marker gene list for `CellAssign`, as well as setting up `scpca-nf` to run with these custom references.
To perform cell type annotation with your own references, you will need to follow these steps:

1. [Create the custom reference(s)](#create-custom-cell-type-references), for use with `SingleR` and/or `CellAssign`
2. [Create a cell type reference metadata file](#create-a-cell-type-reference-metadata-file)
3. [Run the workflow with your custom references](#run-the-workflow-with-custom-references)

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**

- [Create custom cell type references](#create-custom-cell-type-references)
  - [Creating a `SingleR` model from a custom reference dataset](#creating-a-singler-model-from-a-custom-reference-dataset)
    - [Identify an appropriate reference dataset](#identify-an-appropriate-reference-dataset)
    - [Train the `SingleR` model](#train-the-singler-model)
    - [Add additional fields to the model object](#add-additional-fields-to-the-model-object)
      - [Special considerations when using ontology labels](#special-considerations-when-using-ontology-labels)
    - [Export the model](#export-the-model)
  - [Creating a custom `CellAssign` reference from a marker-gene list](#creating-a-custom-cellassign-reference-from-a-marker-gene-list)
- [Create a cell type reference metadata file](#create-a-cell-type-reference-metadata-file)
- [Run the workflow with custom references](#run-the-workflow-with-custom-references)
  - [Special considerations for repeating cell type annotation](#special-considerations-for-repeating-cell-type-annotation)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->



## Create custom cell type references


### Creating a `SingleR` model from a custom reference dataset

Your `SingleR` reference dataset should be a well-curated gene expression dataset that contains samples, or cells in the case of single-cell expression data, with known cell type annotations.
[As described in the `SingleR` book](https://bioconductor.org/books/release/SingleRBook/), `SingleR` uses these reference expression values and associated labels to identify cells in your dataset with similar expression patterns.

The `scpca-nf` workflow requires a _pre-trained `SingleR` model_ as built from your selected reference dataset.
Therefore, you will need to install the `SingleR` package to prepare a reference model:

```
# install BiocManager, as needed, to install Bioconductor packages
# install.packages("BiocManager")

# install SingleR
BiocManager::install("SinglerR")
```

Trained model files should be saved as RDS files and named `<singler_reference_name>_model.rds`, where `<singler_reference_name>` is a string of your choosing.
Create this file by taking the following steps:

1. [Identify an appropriate reference dataset](#identify-an-appropriate-reference-dataset)
2. [Train your model](#train-the-singler-model) based on your reference dataset with [`SingleR::trainSingleR()`](https://rdrr.io/github/LTLA/SingleR/man/trainSingleR.html)
3. [Add additional fields](#add-additional-fields-to-the-model-object) to the established model object
4. [Export the trained model](#export-the-model) to a file named `<singler_reference name>_model.rds`

#### Identify an appropriate reference dataset

Either bulk gene expression (microarray or RNA-seq) or single-cell datasets with annotations can be used as references.
To the best of your ability, you should identify and use a reference dataset that contains cell types that you expect will also be present in your data.
The reference dataset should be from the same organism as the sample you are analyzing.

Your reference dataset should be formatted as either a [`SummarizedExperiment`](https://rdrr.io/bioc/SummarizedExperiment/man/SummarizedExperiment-class.html) or [`SingleCellExperiment`](https://rdrr.io/bioc/SingleCellExperiment/man/SingleCellExperiment.html) object, where each column represents a sample or cell, and each row represents a gene.
Rows should be named with Ensembl gene IDs.
For more information on creating these types of objects, please see the following resources:

- [Constructing a `SummarizedExperiment` object](https://bioconductor.org/packages/release/bioc/vignettes/SummarizedExperiment/inst/doc/SummarizedExperiment.html#constructing-a-summarizedexperiment)
- An introduction to [the `SingleCellExperiment` class](https://bioconductor.org/books/3.18/OSCA.intro/the-singlecellexperiment-class.html#storing-primary-experimental-data)
- [Converting `Seurat` objects](https://satijalab.org/seurat/reference/as.singlecellexperiment) or [`AnnData` objects](https://bioconductor.org/packages/3.19/bioc/html/zellkonverter.html) to `SingleCellExperiment` objects

Some resources you may find helpful for obtaining a suitable reference dataset are the Bioconductor packages [`celldex`](https://bioconductor.org/packages/release/data/experiment/vignettes/celldex/inst/doc/userguide.html) and [`scRNA-seq`](https://bioconductor.org/packages/release/data/experiment/html/scRNAseq.html), and the [`Azimuth`](https://azimuth.hubmapconsortium.org/references/) database.
Note that the [`SingleR` reference datasets which `scpca-nf` has pre-compiled for use](./external-instructions.md#singler-references) are a selected subset of references from the `celldex` package.

#### Train the `SingleR` model

To train your model, use the following code, where

- `ref` is your `SummarizedExperiment` or `SingleCellExperiment` reference dataset object with columns as samples or cells and rows as genes.
- `labels` is a vector providing cell type labels for the sample or cells in the reference dataset.
We encourage using [Cell Ontology IDs for cell type labels](https://www.ebi.ac.uk/ols/ontologies/cl) when performing annotation.
Ontology IDs provide a standardized terminology for cell types and also allow for more complex relationships between cell types.
[See here for additional considerations](#special-considerations-when-using-ontology-labels) if you chose to provide ontology ids.
- `restrict` should be a vector of Ensembl gene IDs corresponding to the genes in the _mapping reference_ used by `scpca-nf`.
You can optionally use a subset of those genes, if you want the model to be focused on a particular set of genes.

```
singler_model <- SingleR::trainSingleR(
  ref = <your reference se or sce>,
  labels = <vector of labels for each sample or cell>,
  restrict = <optionally, vector of gene ids that are in the mapping reference>
)
```

#### Add additional fields to the model object

Once the model object is built, you can optionally add two fields to the model object:

- `reference_label`, a string with information about the types of labels you specified
- `reference_name`, a string with your reference name

If you choose to add this information, use the following code:

```
singler_model$reference_label <- "<name describing your labels>"
singler_model$reference_name <- "<name of your singler model>"
```

If you do not provide this information, the final `SingleCellExperiment` object produced by `scpca-nf` will have `NULL` values for the fields `metadata(sce)$singler_reference_label` and `metadata(sce)$singler_reference_name`, respectively, but no calculations will be affected.


##### Special considerations when using ontology labels

If you chose to use ontology labels, you may wish to set the `reference_label` field to the string `"label.ont"` as shown below.

```
singler_model$reference_label <- "label.ont"
```

If you use `"label.ont"`, you must also include a data frame named `celltype_ontology_df` in the `SingleR` model object, which contains two columns: `ontology_id` and `ontology_cell_names`.
The `ontology_id` column should contain ontology id values corresponding to labels in your `labels` vector, and the `ontology_cell_names` column should contain the human-readable cell type names associated with each ontology ID (e.g., `"B cell"`).
The [`ontoProc` package](https://bioconductor.org/packages/release/bioc/html/ontoProc.html) from `Bioconductor` may be useful for compiling this information.

```
singler_model$celltype_ontology_df <- data.frame(
  ontology_id = <vector of ontology ids corresponding to annotation labels>,
  ontology_cell_names = <vector of matching cell type names>
)
```

#### Export the model

Finally, export the trained model to a file named `<singler_reference_name>_model.rds`:
```
model_file_name <- "<singler_reference_name>_model.rds"
readr::write_rds(singler_model, model_file_name)
```



### Creating a custom `CellAssign` reference from a marker-gene list

The `CellAssign` reference file is created by converting a list of marker genes for a set of cell type labels into a binary matrix with values of `0` and `1`.
This matrix should have all possible cell types as the columns and all possible genes, represented as Ensembl gene IDs, as the rows.
Values of `1` indicate that the given gene is a marker gene for the given cell type, and values of `0` indicate that the gene is not a marker gene for the cell type.
When compiling this information, be sure to use Ensembl gene IDs from the same organism as the sample you are analyzing!

Some resources that you might find helpful for compiling marker gene lists include [PanglaoDB](https://panglaodb.se/) (note that the [pre-compiled `CellAssign` references in `scpca-nf`] were obtained from `PanglaoDB` marker gene lists), [`MSigDB`](https://www.gsea-msigdb.org/gsea/msigdb/genesets.jsp?collection=C8), [`CellMarker`](http://bio-bigdata.hrbmu.edu.cn/CellMarker/), and [`singleCellBase`](http://cloud.capitalbiotech.com/SingleCellBase/).
Note that the [`CellAssign` reference datasets which `scpca-nf` has pre-compiled for use](./external-instructions.md/#cellassign-references) were derived from marker gene sets in `PanglaoDB`.

`CellAssign` reference files should be saved as TSV files and named `<cellassign_reference_name>.tsv`, where `<cellassign_reference_name>` is a string of your choosing.

This TSV file should have the following columns and values:

1. `ensembl_id`, which contains Ensembl gene IDs.
2. One column for each cell type label to use for annotation.
All values in this column should be `1` or `0`, indicating whether that row's Ensembl gene id is a marker gene for that cell type (`1`) or not (`0`)
3. `other`, which contains all values set to `0`.
If `CellAssign` cannot confidently annotate a cell with one of the provided cell type labels, then it will assign the label of `"other"`, indicating that no marker genes were specifically expressed in that cell.

Once this file is created, save it as `<cellassign_reference_name>.tsv`.


## Create a cell type reference metadata file

When performing cell type annotation, `scpca-nf` requires a TSV metadata file with information about cell type references, [as described in these instructions](./external-instructions.md#preparing-the-project-cell-type-metadata-file).

All libraries with the same `scpca_project_id` will use the same reference dataset for each of `SingleR` and `CellAssign`, respectively.
The project cell type metadata file should contain these five columns with the following information:

| column_id             | contents |
| --------------------- | -------- |
| `scpca_project_id`    | Project ID matching values in the run metadata file |
| `singler_ref_name`    | The name of your `SingleR` model which should match the `<singler_reference_name>` string you used when creating your `SingleR` reference RDS file. Use `NA` to skip `SingleR` annotation |
| `singler_ref_file`    | The file name of your `SingleR` model, as `<singler_reference_name>_model.rds`. Use `NA` to skip `SingleR` annotation |
| `cellassign_ref_name` |The name of your `CellAssign` model which should match the `<cellassign_reference_name>` string you used when creating your `CellAssign` reference TSV file. Use `NA` to skip `CellAssign` annotation |
| `cellassign_ref_file` | The file name of your `CellAssign` model, as `<cellassign_reference_name>.tsv`. Use `NA` to skip `CellAssign` annotation |

## Run the workflow with custom references

There are several parameter values you will have to set either in your [configuration file](./external-instructions.md#configuration-files) or at the command line when running `scpca-nf`:


- `params.singler_models_dir`: The path/uri to the directory that contains your `<singler_reference_name>_model.rds`.
This value should _not include_ the file name itself.
For example, if the path to your file is `/path/to/reference/<singler_reference_name>_model.rds`, only provide the value `/path/to/reference`.
You can skip this parameter if you are not performing cell type annotation with `SingleR`.
- `params.cellassign_reference_dir`: The path/uri to the directory that contains your `<cellassign_reference_name>.tsv`.
This value should _not include_ the file name itself.
For example, if the path to your file is `/path/to/reference/<cellassign_reference_name>.tsv`, only provide the value `/path/to/reference`.
You can skip this parameter if you are not performing cell type annotation with `CellAssign`.
- `params.project_celltype_metafile`: The path/uri to your cell type metadata file, including the file name.


When running the workflow, be sure to use the `--perform_celltyping` flag to perform cell type annotation.

You can run the workflow using the below code.
The example below shows all parameters specified at the command line, but if specifying them in the configuration file you would be able to remove those lines.

```bash
nextflow run AlexsLemonade/scpca-nf \
  --singler_models_dir <directory/containing/singler/model/file> \
  --cellassign_reference/dir <directory/containing/cellassign/reference/file> \
  --project_celltype_metafile <path/to/celltype metafile> \
  --perform_celltyping
```

### Special considerations for repeating cell type annotation


As described [in these instructions](./external-instructions.md#repeating-cell-type-annotation), cell type annotation is, by default, not repeated if results already exist.
`scpca-nf` determines if results already exist based on two criteria: i) the _file name_ of the provided reference file has not changed, and ii) [intermediate checkpoint files](./external-instructions.md#repeating-cell-type-annotation) from previous runs of cell typing already exist.
For example, if you re-run the workflow with a renamed but unchanged reference file, cell type annotation will be repeated.

To force cell type annotation to repeat, for example if you retain the same reference file names but changed their contents, use the `--repeat_celltyping` flag as well as the `--perform_celltyping` flag at runtime:

```bash
nextflow run AlexsLemonade/scpca-nf \
  --singler_models_dir <directory/containing/singler/model/file> \
  --cellassign_reference/dir <directory/containing/cellassign/reference/file> \
  --project_celltype_metafile <path/to/celltype metafile> \
  --perform_celltyping \
  --repeat_celltyping
```
