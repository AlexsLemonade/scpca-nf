#!/bin/bash

# create a new, detached, session named "nextflow" and run the run_nextflow script in it
tmux new-session -d -s "nextflow" "/opt/nextflow/scripts/run_nextflow.sh"
