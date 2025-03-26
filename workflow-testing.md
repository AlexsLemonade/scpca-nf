# Testing scpca-nf outputs

As `scpca-nf` is updated, we want to be sure that we are aware of any changes in the output from the pipeline before we deploy updated results to the portal.
This document outlines the tests that we will run and the key metrics that we will monitor to ensure that the pipeline is functioning as expected.

## Data format checks

Data formatting checks will be conducted within the `scpca-nf` workflow to ensure that all output files are in the expected format.
Details to come.

## Output metrics

The table below describes the key metrics that we will track for the workflow outputs, and which files are required to calculate those metrics.
These metrics should be calculated, if required, and compiled at the end of the workflow just before the publication step, and stored in a `_metrics.json` file.

| Metric                       | Description                                                       | Required files                                                |
| ---------------------------- | ----------------------------------------------------------------- | ------------------------------------------------------------- |
| `unfiltered_cells`           | total number of cells in unfiltered SCE object                    | `{library_id}_metadata.json` or `{library_id}_unfiltered.rds` |
| `filtered_cells`             | total number of cells in filtered SCE object                      | `{library_id}_metadata.json` or `{library_id}_filtered.rds`   |
| `processed_cells`            | total number of cells in processed SCE object                     | `{library_id}_metadata.json` or `{library_id}_processed.rds`  |
| `total_reads`                | total number of reads from FASTQ files                            | `{library_id}_metadata.json`                                  |
| `mapped_reads`               | total number of reads mapped by alevin-fry                        | `{library_id}_metadata.json`                                  |
| `droplet_filtering_method`   | method used to filter droplets                                    | `{library_id}_metadata.json`                                  |
| `normalization_method`       | method used to normalize counts                                   | `{library_id}_metadata.json` or `{library_id}_processed.rds`  |
| `cell_filtering_method`      | method used to filter cells                                       | `{library_id}_metadata.json` or `{library_id}_processed.rds`  |
| `miqc_pass_count`            | count of cells that pass miQC filtering                           | `{library_id}_filtered.rds`                                   |
| `unfiltered_total_counts`    | sum of counts matrix in unfiltered SCE object                     | `{library_id}_unfiltered.rds`                                 |
| `unfiltered_total_spliced`   | sum of spliced counts matrix in unfiltered SCE object             | `{library_id}_unfiltered.rds`                                 |
| `filtered_total_counts`      | sum of counts matrix in filtered SCE object                       | `{library_id}_filtered.rds`                                   |
| `filtered_total_spliced`     | sum of spliced counts matrix in filtered SCE object               | `{library_id}_filtered.rds`                                   |
| `processed_total_counts`     | sum of counts matrix in processed SCE object                      | `{library_id}_processed.rds`                                  |
| `processed_total_spliced`    | sum of spliced counts matrix in processed SCE object              | `{library_id}_processed.rds`                                  |
| `processed_total_logcounts`  | sum of logcounts matrix in processed SCE object                   | `{library_id}_processed.rds`                                  |
| `unfiltered_expressed_genes` | count of genes with total count > 0 in unfiltered SCE object      | `{library_id}_unfiltered.rds`                                 |
| `filtered_expressed_genes`   | count of genes with total count > 0 in filtered SCE object        | `{library_id}_filtered.rds`                                   |
| `processed_expressed_genes`  | count of genes with total count > 0 in processed SCE object       | `{library_id}_processed.rds`                                  |
| `hv_genes`                   | the array of the highly variable genes                            | `{library_id}_processed.rds`                                  |
| `cluster_algorithm`          | the algorithm used to cluster cells                               | `{library_id}_processed.rds`                                  |
| `cluster_sizes`              | an array with the number of cells in each cluster                 | `{library_id}_processed.rds`                                  |
| `singler_reference`          | The reference used for SingleR                                    | `{library_id}_processed.rds`                                  |
| `singler_celltypes`          | a dictionary with cell counts for each SingleR-based cell type    | `{library_id}_processed.rds`                                  |
| `cellasssign_reference`      | the reference used for CellAssign                                 | `{library_id}_processed.rds`                                  |
| `cellassign_celltypes`       | a dictionary with cell counts for each CellAssign-based cell type | `{library_id}_processed.rds`                                  |
| `unfiltered_altexp_total`    | total number of reads from unfiltered altExps (dict of values)    | `{library_id}_unfiltered.rds`                                 |
| `filtered_altexp_total`      | total number of reads from filtered altExps (dict of values)      | `{library_id}_filtered.rds`                                   |
| `processed_altexp_total`     | total number of reads from processed altExps (dict of values)     | `{library_id}_processed.rds`                                  |
| `adt_scpca_filter_count`     | count of cells that pass ADT filtering (labeled as `Keep`)        | `{library_id}_filtered.rds`                                   |

### Metric comparisons

After a workflow run is complete, we can use the `_metrics.json` files to compare the metrics from the current run to the metrics from a previous run.
We will write a script to perform these comparisons and generate a report that highlights any significant changes in the metrics.
This will allow us to identify any changes in the output metrics that may indicate unexpected effects of changes to the pipeline.
