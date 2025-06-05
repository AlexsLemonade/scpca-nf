#!/usr/bin/env python3

# This script creates the required configuration file for running cellranger multi

from pathlib import Path
import argparse
import textwrap
import re

parser = argparse.ArgumentParser()
parser.add_argument(
    "--config", required=True, type=Path, help="Path to output config CSV file"
)
parser.add_argument(
    "--transcriptome_reference",
    required=True,
    type=Path,
    help="Path to transcriptome reference",
)
parser.add_argument(
    "--probe_set_reference",
    required=True,
    type=Path,
    help="Path to probe set reference",
)
parser.add_argument(
    "--fastq_dir",
    required=True,
    type=Path,
    help="Path to directory containing input FASTQ files",
)

args = parser.parse_args()

# check that path to reference files exist
if not args.transcriptome_reference.exists():
    raise FileNotFoundError(
        f"Transcriptome reference not found: {args.transcriptome_reference}"
    )

if not args.probe_set_reference.exists():
    raise FileNotFoundError(
        f"Probe set reference not found: {args.probe_set_reference}"
    )

# check that path to FASTQ directory exists
if not args.fastq_dir.exists():
    raise FileNotFoundError(f"FASTQ directory not found: {args.fastq_dir}")

# build config file content
config_content = textwrap.dedent(
    f"""
    [gene-expression]
    reference,{args.transcriptome_reference.resolve()}
    probe-set,{args.probe_set_reference.resolve()}
    create-bam,false
    
    [libraries]
    fastq_id,fastqs,feature_types
    """.lstrip()  # remove leading newline
)

# regex for extracting ID from the FASTQ files
pattern = re.compile(r"^(.+)_S.+_L.+_[RI].+\.fastq\.gz$")

# Get sample IDs needed for cellranger from fastq files
fastq_path = args.fastq_dir.resolve()
fastq_ids = set()
for fastq_file in fastq_path.glob("*.fastq.gz"):
    match = pattern.match(fastq_file.name)
    if match:
        fastq_ids.add(match.group(1))

# make sure sample ids were found
if not fastq_ids:
    raise ValueError(
        f"{fastq_path} does not contain FASTQ files that match the expected naming convention."
    )

# add fastq_ids to config content
for fastq_id in fastq_ids:
    config_content += f"{fastq_id},{fastq_path},Gene_Expression\n"

# save config content to file
args.config.write_text(config_content)
