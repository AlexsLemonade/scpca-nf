#!/bin/bash
set -euo pipefail
# This script is used to build and/or push the Cavatica app

# Run from the script file location
cd "$(dirname "${BASH_SOURCE[0]}")"

# cavatica app id: <username>/<project>/<app>
app_id="jashapiro/scpca-nf-test/scpca-nf"
doc_file="sb_doc.md"
exclude_file="sb_exclude.txt"

# get the action from the command line, default is build
ACTION=${1:-"build"}

if [ "$ACTION" == "build" ]; then
    # default to multi-instance mode
    mode_arg="--execution-mode multi-instance"

    # use previous schema if it exists (moving current to backup)
    if [ -f sb_nextflow_schema.yaml ]; then
        mv sb_nextflow_schema.yaml sb_nextflow_schema.yaml.bak
        #replace execution-mode with previous schema
        mode_arg="--sb-schema sb_nextflow_schema.yaml.bak"
    fi

    # Build sb_nextflow_schema.yaml
    pixi run -e cavatica sbpack_nf \
      --workflow-path .. \
      --appid ${app_id} \
      --sb-doc ${doc_file} \
      ${mode_arg} \
      --dump-sb-app
    mv ../sb_nextflow_schema.yaml .
elif [ "$ACTION" == "push" ]; then
    # Check if the sb_nextflow_schema.yaml file exists
    if [ ! -f sb_nextflow_schema.yaml ]; then
        echo "sb_nextflow_schema.yaml file not found. Please run the build action first."
        exit 1
    fi
    # Push the Docker image to the registry
    pixi run -e cavatica sbpack_nf \
      --profile . \
      --workflow-path .. \
      --sb-schema sb_nextflow_schema.yaml \
      --appid ${app_id} \
      --exclude $(< $exclude_file)
else
    echo "Invalid action: $ACTION"
    exit 1
fi
