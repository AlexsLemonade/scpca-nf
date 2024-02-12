# Overview

The `scpca-nf` workflow includes the option to perform cell type annotation, [as described in these instructions](./external-instructions.md#cell-type-annotation), with the reference-based method [`SingleR`](https://bioconductor.org/packages/release/bioc/html/SingleR.html) and/or the marker-gene based method [`CellAssign`](https://github.com/Irrationone/cellassign).

Several [reference datasets](./external-instructions.md#singler-references) and [marker gene lists](./external-instructions.md#cellassign-references) have been pre-compiled for use.

This document provides instructions for creating and using your own reference data set for `SingleR` annotation or marker gene list for `CellAssign`, as well as setting up `scpca-nf` to run with these custom references.
To perform cell type annotation with your own references, you will need to follow these steps:

1. Create the custom reference(s), for use with `SingleR` and/or `CellAssign`
2. Create a cell type reference metadata file
3. Run the workflow with your custom references

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**

- [Create custom cell type references](#create-custom-cell-type-references)
  - [Using a custom `SingleR` reference dataset](#using-a-custom-singler-reference-dataset)
    - [Training the model](#training-the-model)
    - [Special considerations when using ontology labels](#special-considerations-when-using-ontology-labels)
  - [Creating a custom `CellAssign` marker-gene list](#creating-a-custom-cellassign-marker-gene-list)
- [Create a cell type reference metadata file](#create-a-cell-type-reference-metadata-file)
- [Run the workflow with custom references](#run-the-workflow-with-custom-references)
  - [Special considerations for repeating cell type annotation](#special-considerations-for-repeating-cell-type-annotation)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->



## Create custom cell type references


### Using a custom `SingleR` reference dataset

Your `SingleR` reference dataset should be a well-curated gene expression dataset that contains samples, or cells in the case of single-cell expression data, with known cell type annotations.
[As described in the `SingleR` book](https://bioconductor.org/books/release/SingleRBook/), `SingleR` uses these reference expression values and associated labels to identify cells in your dataset with similar expression patterns.
Either bulk gene expression (microarray or RNA-seq) or single-cell datasets with annotations can be used as references.
To the best of your ability, you should identify and use a reference dataset that contains cell types that you expect will also be present in your data.
The reference dataset should also be the same organism as the sample you are analyzing to ensure matching gene ids.

The `scpca-nf` workflow required a _pre-trained `SingleR` model_ as built from your selected reference dataset.
These trained model files should be saved as RDS files and named `<singler_reference_name>_model.rds`, where `<singler_reference_name>` is a string of your choosing.
Create this file by taking the following steps:

1. Identify an appropriate reference dataset to use.
This dataset should be either a [`SummarizedExperiment`](https://rdrr.io/bioc/SummarizedExperiment/man/SummarizedExperiment-class.html) or [`SingleCellExperiment`](https://rdrr.io/bioc/SingleCellExperiment/man/SingleCellExperiment.html) object, where each column represents a sample or cell, and each row represents a gene.
Rows should be named with Ensembl gene ids.
  - Some resources you can use to find a suitable reference dataset are the [`celldex` package](https://bioconductor.org/packages/release/data/experiment/vignettes/celldex/inst/doc/userguide.html) or the [`scRNA-seq` package](https://bioconductor.org/packages/release/data/experiment/html/scRNAseq.html), both from Bioconductor.
Note that the [`SingleR` reference datasets which `scpca-nf` has pre-compiled for use](./external-instructions.md#singler-references) are a selected subset of references from the `celldex` package.
2. Use [`SingleR::traingSingleR()`](https://rdrr.io/github/LTLA/SingleR/man/trainSingleR.html) to train your model based on your chosen reference dataset.
3. Add additional fields to the established model object, as needed (see instructions below).
4. Export the trained model to a file `<singler_reference name>_model.rds`.


#### Training the model

To train your model, use the following code, where

- `ref` is your `SummarizedExperiment` or `SingleCellExperiment` reference dataset object with columns as samples or cells and rows as genes.
- `labels` is a vector providing cell type labels for the sample or cells in the reference dataset.
We strongly encourage the use of ontology cell type labels when performing annotation; [see here for additional considerations](#special-considerations-for-using-ontology-labels) if you chose to provide ontology ids.
- `restrict` can optionally be used to only consider genes which are present in the _mapping reference_ used by `scpca-nf`.
This should be a vector of Ensembl gene ids.

```
singler_model <- SingleR::trainSingleR(
  ref = <your reference se or sce>,
  labels = <vector of labels for each sample or cell>,
  restrict = <optionally, vector of gene ids that are in the mapping reference>
)
```

Once the model object is built, you can optionally add the field `reference_label` with information about the types of labels you specified.
If you do not provide this information, the final `SingleCellExperiment` object produced by `scpca-nf` will have a value `NULL` for the field `metadata(sce)$singler_reference_label`, but no calculations will be affected.

If you choose to add this label, use the following code:

```
singler_model$reference_label <- "<name describing my labels>"
```

Finally, export the trained model to a file named `<singler_reference_name>_model.rds`:

```
model_file_name <- "<singler_reference_name>_model.rds"
readr::write_rds(singler_model, model_file_name)
```


#### Special considerations when using ontology labels

If you chose to use ontology labels, you may wish to set the `reference_label` field to the string `"label.ont"`.
Using this specific string will trigger the `scpca-nf` workflow to additionally expect a dataframe with two columns, `ontology_id` and `ontology_cell_names`, to be present in the `SingleR` model object in the field `celltype_ontology_df`.
The `ontology_id` column should contain ontology id values corresponding to labels in your `labels` vector, and the `ontology_cell_names` column should contain the human-readable cell type names associated with each ontology ID (e.g., `"B cell"`).
The [`oncoProc` package](https://bioconductor.org/packages/release/bioc/html/ontoProc.html) from `Bioconductor` may be useful for compiling this information.

```
singler_model$reference_label <- "label.ont"
singler_model$celltype_ontology_df <- data.frame(
  ontology_id = <vector of ontology ids corresponding to annotation labels>,
  ontology_cell_names = <vector of matching cell type names>
)
```

### Creating a custom `CellAssign` marker-gene list

`CellAssign` reference files should be saved as TSV files and named `<cellassign_reference_name>.tsv`, where `<cellassign_reference_name>` is a string of your choosing.

This TSV file should have the following columns and values:

1. `ensembl_id`, which contain Ensembl gene ids.
2. One column for each cell type label to use for annotation.
All values in this column should be `1` or `0`, indicating whether that row's Ensembl gene id is a marker gene for that cell type (`1`) or not (`0`)
3. Optionally, you may wish to a include a column called `"other"` with all values set to `0`.
If `CellAssign` cannot confidently annotate a cell with one of the provided cell type labels, then it will assign the label of `"other"`, indicating that no marker genes were specifically expressed in that cell.

Once this file is created, save it as `<cellassign_reference_name>.tsv`.


## Create a cell type reference metadata file

When performing cell type annotation, `scpca-nf` requires a TSV metadata file with information about cell type references, [as described in these instructions](./external-instructions.md#preparing-the-project-cell-type-metadata-file).

All libraries a given `scpca_project_id` will use the same reference dataset for each of `SingleR` and `CellAssign`, respectively.
The project cell type metadata file should contain these five columns with the following information:

| column_id             | contents |
| --------------------- | -------- |
| `scpca_project_id`    | Project id matching values in the run metadata file |
| `singler_ref_name`    | The name of your `SingleR` model which should match the `<singler_reference_name>` string you used when creating your `SingleR` reference RDS file. Use `NA` to skip `CellAssign` annotation |
| `singler_ref_file`    | The file name of your `SingleR` model, as `<singler_reference_name>_model.rds`. Use `NA` to skip `CellAssign` annotation |
| `cellassign_ref_name` |The name of your `CellAssign` model which should match the `<cellassign_ref_name>` string you used when creating your `CellAssign` reference TSV file. Use `NA` to skip `CellAssign` annotation |
| `cellassign_ref_file` | The file name of your `CellAssign` model, as `<cellassign_reference_name>.tsv`. Use `NA` to skip `CellAssign` annotation |

## Run the workflow with custom references

There are several parameter values you will have to set either in your [configuration file](./external-instructions.md#configuration-files) or at the command line when running `scpca-nf`:


- `params.singler_reference_label`: A name describing the label used during `SingleR` model training.
If you do not have a particular name for your labels, you can simply provide `""`, _but you MUST provide a non-default value for this parameter_ unless you are not performing cell type annotation with `SingleR`.
- `params.singler_models_dir`: The path/uri to the directory that contains your `<singler_reference_name>_model.rds`.
This value should _not include_ the file name itself.
For example, if the path to your file is `/path/to/reference/<singler_reference_name>_model.rds`, only provide the value `/path/to/reference`.
You can skip this parameter if you are not performing cell type annotation with `SingleR`.
- `params.cellassign_reference_dir`: The path/uri to the directory that contains your `<cellassign_reference_name>.tsv`.
This value should _not include_ the file name itself.
For example, if the path to your file is `/path/to/reference/<cellassign_reference_name>.tsv`, only provide the value `/path/to/reference`.
You can skip this parameter if you are not performing cell type annotation with `CellAssign`.
- `params.project_celltype_metafile`: The path/uri to your cell type metadata file, including the the file name.


When running the workflow, be sure to use the `--perform_celltyping` flag to perform cell type annotation.

You can run the workflow as (showing all parameters specified at the command line, for example):

```bash
nextflow run AlexsLemonade/scpca-nf \
  --singler_label_name <optionally, the label name> \
  --singler_models_dir <directory/containing/singler/model/file> \
  --cellassign_reference/dir <directory/containing/cellassign/reference/file> \
  --project_celltype_metafile <path/to/celltype metafile> \
  --perform_celltyping
```

### Special considerations for repeating cell type annotation


As described [in these instructions](./external-instructions.md#repeating-cell-type-annotation), cell type annotation is, by default, not repeated if results already exist.
`scpca-nf` determines if results already exist based on the _file name_ of the provided reference file.
Therefore, if you re-run the workflow with a renamed, but unchanged, reference file, be aware that cell type annotation will be repeated.
To force cell type annotation to repeat, for example if you retain the same reference file names but changed their contents, use the `--repeat_celltyping` flag as well as the `--perform_celltyping` flag at runtime.