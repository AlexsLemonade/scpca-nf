
## Comparing results between `scpca-nf` runs

To facilitate comparisons between the production version of `scpca-nf` results and the results from a staging (or other) run, we have written a script and notebook to compare metrics between different runs, producing an HTML report.

The script can be run for a single project using a command like the following, assuming that AWS credentials are available in the environment:

```
Rscript compare-metrics.R \
  --project_id SCPCP000001 \
  --output_file "SCPCP000001_metrics_comparison.html"
```

Multiple projects can be included using a comma-separated list for the `--project_id` argument.

The default comparison will be between the metrics files available in the production and staging directories, treating the production versions as the reference run.
The files are expected to be found at `s3://ccdl-scpca-results-prod-997241705947-us-east-1/results` and `s3://ccdl-scpca-results-staging-997241705947-us-east-1/results`, by default, but these can be changed using the `--ref_s3` and `--comp_s3` arguments.
