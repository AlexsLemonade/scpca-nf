# Testing scpca-nf outputs

As `scpca-nf` is updated, we want to be sure that we are aware of any changes in the output from the pipeline before we deploy updated results to the portal.
This documents outlines the tests that we will run and the key metrics that we will monitor to ensure that the pipeline is functioning as expected.

## Data format checks

## Output metrics

The table below describes the key metrics that we will track for the workflow outputs, and which files are required to calculate those metrics.
These metrics should be calculated, if required, at the end of the workflow just before the publication step, and stored in a JSON file. At the moment, I am not sure if we just add these to `_metadata.json` or create a separate file.

| Metric             | Description                                    | Required files              |
| ------------------ | ---------------------------------------------- | --------------------------- |
| `unfiltered_cells` | total number of cells in unfiltered SCE object | `{library_id}_metadata.json` or `{library_id}_unfiltered.rds` |
| `filtered_cells`   | total number of cells in filtered SCE object   | `{library_id}_metadata.json` or `{library_id}_filtered.rds` |
| `processed_cells`  | total number of cells in processed SCE object  | `{library_id}_metadata.json` or `{library_id}_processed.rds` |
| `total_reads`      | total number of reads from FASTQ files | `{library_id}_metadata.json` |
| `mapped_reads`      | total number of reads mapped by alevin-fry | `{library_id}_metadata.json` |
| `cell_filtering_method` | method used to filter cells | `{library_id}_metadata.json` or `{library_id}_processed.rds` |
| `miqc_pass_count` | count of cells that pass miQC filtering | `{library_id}_filtered.rds` |
| `unfiltered_total_expression` | total expression in unfiltered SCE object (sum of counts matrix) | `{library_id}_unfiltered.rds` |
| `filtered_total_expression` | total expression in filtered SCE object (sum of counts matrix) | `{library_id}_filtered.rds` |
| `processed_total_expression` | total expression in processed SCE object (sum of counts matrix) | `{library_id}_processed.rds` |
| `unfiltered_expressed_genes` | count of genes with total expression > 0 in unfiltered SCE object | `{library_id}_unfiltered.rds` |
| `filtered_expressed_genes` | count of genes with total expression > 0 in filtered SCE object | `{library_id}_filtered.rds` |
| `processed_expressed_genes` | count of genes with total expression > 0 in processed SCE object | `{library_id}_processed.rds` |
| `hv_genes` | the array of the highly variable genes  | `{library_id}_processed.rds` |
| `cluster_sizes` | an array with the number of cells in each cluster | `{library_id}_processed.rds` |
| `singler_celltypes` | a dictionary with the number of cells in each SingleR-based cell type  | `{library_id}_processed.rds` |
| `cellassign_celltypes` | a dictionary with the number of cells in each CellAssign-based cell type  | `{library_id}_processed.rds` |
| `unfiltered_altexp_total` | total number of reads from altExp (ADT or cellhash) | `{library_id}_unfiltered.rds` |
| `filtered_altexp_total` | total number of reads from altExp (ADT or cellhash) | `{library_id}_filtered.rds` |
| `processed_altexp_total` | total number of reads from altExp (ADT or cellhash) | `{library_id}_processed.rds` |
