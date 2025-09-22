This directory contains information and metadata about references used in the `scpca-nf` workflow.

If you are a Data Lab member looking for more information about maintaining or adding additional references, please refer to [`internal-instructions.md`](../internal-instructions.md).

### Transcriptome reference files

These files are used by `build-index.nf` to generate transcriptome reference files.

- `ref-metadata.tsv` provides transcriptome reference metadata for organisms considered in the `scpca-nf` workflow.
This file is used by `build-index.nf` to generate transcriptome references.
- `scpca-refs.json` contains the relative paths of the default reference files for mapping and quantification provided by the Data Lab, and is used for the `scpca-nf` workflow configuration parameter `ref_json`.
  The base location for the reference files is given by the configuration parameter `ref_rootdir`.
  This file is produced by `scripts/create-reference-json.R`.


### Cell type reference files

These files are used by `build-celltype-ref.nf` to generate cell type annotation reference files.

- `celltype-reference-metadata.tsv` provides information about the name, source, and contents of cell type annotation references established by the Data Lab for use with `SingleR` and `CellAssign`.
- `PanglaoDB_markers_<date>.tsv` is the full set of marker genes exported from [PanglaoDB](https://panglaodb.se), versioned at the given date.

These files are copies of reference files that were created as part of `OpenScPCA`.
See the [`cell-type-consensus` module as part of the `OpenScPCA-analysis` repository](https://github.com/AlexsLemonade/OpenScPCA-analysis/tree/v0.2.3/analyses/cell-type-consensus) to learn more about how these file were originally created.
These files were copied from the `v0.2.3` release. 

- [`panglao-cell-type-ontologies.tsv`](https://github.com/AlexsLemonade/OpenScPCA-analysis/blob/v0.2.3/analyses/cell-type-consensus/references/panglao-cell-type-ontologies.tsv) contains the [Cell Ontology identifier term](https://www.ebi.ac.uk/ols4/ontologies/cl) for all cell types available in the `PanglaoDB` reference used when running `CellAssign`.
The table includes the following columns:

|  |   |
| --- | --- |
| `ontology_id` | [cell type (CL) ontology identifier term](https://www.ebi.ac.uk/ols4/ontologies/cl) |
| `human_readable_value` | Label associated with the cell type ontology term |
| `panglao_cell_type` | Original name for the cell type as set by `PanglaoDB` |

- [`scimilarity-mapped-ontologies.tsv`](https://github.com/AlexsLemonade/OpenScPCA-analysis/blob/v0.2.3/analyses/cell-type-scimilarity/references/scimilarity-mapped-ontologies.tsv) contains the [Cell Ontology identifier term](https://www.ebi.ac.uk/ols4/ontologies/cl) for all cell types available in the `SCimilarity` foundational model, [`model_v1.1`](https://zenodo.org/records/10685499).
The table includes the following columns:

| | | 
| -| - | 
| `scimilarity_celltype_annotation` | The cell type annotation present in the `SCimilarity` model | 
| `cl_annotation` | The human readable name of the cell type found in the Cell Ontology | 
| `scimilarity_celltype_ontology` | The ontology identifier associated with the name of the cell type in the Cell Ontology | 

- [`consensus-cell-type-reference.tsv`](https://github.com/AlexsLemonade/OpenScPCA-analysis/blob/v0.2.3/analyses/cell-type-consensus/references/consensus-cell-type-reference.tsv) contains a table with all cell type combinations between the `PanglaoDB` reference (used with `CellAssign`), the `BlueprintEncodeData` reference (used with `SingleR`), and the `SCimilarity` reference (used with `SCimilarity`) for which a consensus cell type is identified. 

The table includes the following columns: 

|  |   |
| --- | --- |
| `panglao_ontology` | Cell type ontology term for `PanglaoDB` cell type |
| `panglao_annotation` | Name for the `panglao_ontology` term as defined by [Cell Ontology](https://www.ebi.ac.uk/ols4/ontologies/cl) |
| `original_panglao_name` | Original name for the cell type as set by `PanglaoDB` |
| `blueprint_ontology` | Cell type ontology term for `BlueprintEncodeData` cell type |
| `blueprint_annotation_cl` | Name for the `blueprint_ontology` term as defined by [Cell Ontology](https://www.ebi.ac.uk/ols4/ontologies/cl) |
| `scimilarity_ontology` | Cell type ontology term for `SCimilarity` cell type (obtained from `scimilarity-mapped-ontologies.tsv`) |
| `scimilarity_annotation` | Name for the `scimilarity_ontology` term as defined by [Cell Ontology](https://www.ebi.ac.uk/ols4/ontologies/cl) |
| `original_scimilarity_name` | Original name for the cell type as set by `SCimilarity` |
| `cellassign_singler_pair_ontology` | Cell type ontology term for the consensus cell type between `PanglaoDB` (`CellAssign`) and `BlueprintEncodeData` (`SingleR`) |
| `cellassign_singler_pair_annotation` | Name for the consensus cell type between `PanglaoDB` (`CellAssign`) and `BlueprintEncodeData` (`SingleR`) |
| `cellassign_scimilarity_pair_ontology` | Cell type ontology term for the consensus cell type between `PanglaoDB` (`CellAssign`) and `SCimilarity` |
| `cellassign_scimilarity_pair_annotation` | Name for the consensus cell type between `PanglaoDB` (`CellAssign`) and `SCimilarity` |
| `singler_scimilarity_pair_ontology` | Cell type ontology term for the consensus cell type between `SCimilarity` and `BlueprintEncodeData` (`SingleR`) |
| `singler_scimilarity_pair_annotation` | Name for the consensus cell type between `SCimilarity` and `BlueprintEncodeData` (`SingleR`) |
| `consensus_ontology` | Cell type ontology term for consensus cell type |
| `consensus_annotation` | Name for the `consensus_ontology` term as defined by [Cell Ontology](https://www.ebi.ac.uk/ols4/ontologies/cl) |

- [`consensus-validation-groups.tsv`](https://github.com/AlexsLemonade/OpenScPCA-analysis/blob/v0.2.3/analyses/cell-type-consensus/references/consensus-validation-groups.tsv) contains a table of all possible consensus cell type labels and assigns a broader group to use for validating consensus cell type assignments using marker genes.
This table includes the following columns:


|  |   |
| --- | --- |
| `consensus_ontology` | Cell type ontology term for consensus cell type |
| `consensus_annotation` | Human readable name for the consensus cell type |
| `validation_group_ontology` | Cell type ontology term for broader cell type group used for validation |
| `consensus_annotation` | Human readable name for broader cell type group used for validation |

- [`validation-markers.tsv`](https://github.com/AlexsLemonade/OpenScPCA-analysis/blob/v0.2.3/analyses/cell-type-consensus/references/validation-markers.tsv) contains a list of marker genes to aid in validating consensus cell types.
These markers were obtained by grabbing the most frequently observed markers for each cell type included in `consensus-validation-groups.tsv` from the [`CellMarker2.0` list of human marker genes](http://117.50.127.228/CellMarker/CellMarker_download.html).
The top 10 genes (sometimes more if there is a tie in the frequency) are included for each cell type.

The table includes the following columns:

|  |   |
| --- | --- |
| `validation_group_ontology` | Cell type ontology term for broader cell type group used for validation |
| `consensus_annotation` | Human readable name for broader cell type group used for validation |
| `ensembl_gene_id` | Ensembl gene identifier for the marker gene |
| `gene_symbol` | Gene symbol for the marker gene |
| `number_of_tissues` | Total number of tissues that express that marker gene in the specified cell type |
| `celltype_total_tissues` | Total number of tissues that contained the specified cell type |
| `percent_tissues` | Percentage of tissues that express the marker gene in the specified cell type |

### Diagnosis mapping files

These files are used by the main workflow to determine which cells to include in a normal reference for `inferCNV`.

- `broad-diagnosis-map.tsv` maps broad diagnosis groups to individual sample (e.g., submitted) diagnoses in ScPCA.

The table includes the following columns:

|  |   |
| --- | --- |
| `ontology_id` | Ontology term for submitted diagnosis  |
| `human_readable_value` | Human readable name for submitted diagnosis ontology term |
| `submitted_diagnosis` | Submitted diagnosis as recorded in ScPCA objects |
| `diagnosis_group` | Broad diagnosis group |


- `diagnosis-celltype-groups.tsv` maps broad diagnosis groups included in `broad-diagnosis-map.tsv` to consensus cell type validation groups included in `consensus-validation-groups.tsv`.

The table includes the following columns:

|  |   |
| --- | --- |
| `diagnosis_group` | Broad diagnosis group |
| `celltype_groups` | Consensus cell type validation groups to include in an `inferCNV` normal reference for samples of the given broad diagnosis group  |
