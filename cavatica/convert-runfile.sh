#!/bin/bash
set -euo pipefail

project_id="jashapiro/scpca-nf-test"
sample_sheet="example_run_metadata.tsv"
output_file="example_run_metadata-cavatica.tsv"

pixi run -e cavatica sbmanifest \
  --profile . \
  --projectid "$project_id" \
  --sample-sheet "$sample_sheet" \
  --output "$output_file" \
  --columns files_directory feature_barcode_file submitter_cell_types_file

# transform the vs:// paths back to NA if needed (use a temp file to avoid in-place editing)
temp_file=$(mktemp)
sed 's|vs:///Projects/[^/]*/NA|NA|g' "$output_file" > "$temp_file"
mv "$temp_file" "$output_file"
