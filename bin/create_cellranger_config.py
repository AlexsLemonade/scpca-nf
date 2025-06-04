#!/usr/bin/env python3

# This script creates the required configuration file for running cellranger multi

from pathlib import Path
import argparse
import textwrap

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
    "--sample_id",
    required=True,
    type=str,
    help="cellranger sample ID obtained from FASTQ file names",
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
    {args.sample_id},{args.fastq_dir.resolve()},Gene Expression
""".lstrip()
)  # remove leading newline

# save config content to file
args.config.write_text(config_content)
