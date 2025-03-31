This directory contains information and metadata about references used in the `scpca-nf` workflow.

If you are a Data Lab member looking for more information about maintaining or adding additional references, please refer to [`internal-instructions.md`](../internal-instructions.md).

### Transcriptome reference files

These files are used by `build-index.nf` to generate transcriptome reference files.

- `ref-metadata.tsv` provides transcriptome reference metadata for organisms considered in the `scpca-nf` workflow.
This file is used by `build-infex.nf` to generate transcriptome references.
- `scpca-refs.json` contains the relative paths of the default reference files for mapping and quantification provided by the Data Lab, and is used for the `scpca-nf` workflow configuration parameter `ref_json`.
  The base location for the reference files is given by the configuration parameter `ref_rootdir`.
  This file is produced by `scripts/create-reference-json.R`.


### Cell type reference files

These files are used by `build-celltype-ref.nf` to generate cell type annotation reference files.

- `celltype-reference-metadata.tsv` provides information about the name, source, and contents of cell type annotation references established by the Data Lab for use with `SingleR` and `CellAssign`.
- `PanglaoDB_markers_<date>.tsv` is the full set of marker genes exported from [PanglaoDB](https://panglaodb.se), versioned at the given date.

These files are copies of reference files that were created as part of `OpenScPCA`. 
See the [`cell-type-consensus` module as part of the `OpenScPCA-analysis` repository](https://github.com/AlexsLemonade/OpenScPCA-analysis/tree/v0.2.2/analyses/cell-type-consensus) to learn more about how these file were originally created. 

- `panglao-cell-type-ontologies.tsv` contains the [Cell Ontology identifier term](https://www.ebi.ac.uk/ols4/ontologies/cl) for all cell types available in the `PanglaoDB` reference used when running `CellAssign`. 
The table includes the following columns:

|  |   |
| --- | --- |
| `ontology_id` | [cell type (CL) ontology identifier term](https://www.ebi.ac.uk/ols4/ontologies/cl) |
| `human_readable_value` | Label associated with the cell type ontology term |
| `panglao_cell_type` | Original name for the cell type as set by `PanglaoDB` |

- `consensus-cell-type-reference.tsv` contains a table with all cell type combinations between the `PanglaoDB` reference and `BlueprintEncodeData` reference for which a consensus cell type is identified.  
The table includes the following columns: 

|  |   |
| --- | --- |
| `panglao_ontology` | Cell type ontology term for `PanglaoDB` cell type |
| `panglao_annotation` | Name for the `panglao_ontology` term as defined by [Cell Ontology](https://www.ebi.ac.uk/ols4/ontologies/cl) |
| `original_panglao_name` | Original name for the cell type as set by `PanglaoDB` |
| `blueprint_ontology` | Cell type ontology term for `BlueprintEncodeData` cell type |
| `blueprint_annotation_cl` | Name for the `blueprint_ontology` term as defined by [Cell Ontology](https://www.ebi.ac.uk/ols4/ontologies/cl) |
| `consensus_ontology` | Cell type ontology term for consensus cell type |
| `consensus_annotation` | Name for the `consensus_ontology` term as defined by [Cell Ontology](https://www.ebi.ac.uk/ols4/ontologies/cl) |
