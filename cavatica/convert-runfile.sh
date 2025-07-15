#!/bin/bash
set -euo pipefail

pixi run -e cavatica sbmanifest \
  --profile . \
  --projectid "jashapiro/scpca-nf-test" \
  --sample-sheet example_run_metadata.tsv \
  --output example_run_metadata-sb.tsv \
  --columns files_directory \
  --validate
