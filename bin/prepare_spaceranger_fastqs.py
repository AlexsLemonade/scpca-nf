#!/usr/bin/env python3
# prepare_spaceranger_fastqs.py
# Usage: prepare_spaceranger_fastqs.py <fastq_dir> <sample_name> <staged_dir>
#
# Stages allowed FASTQ files for Space Ranger by creating a
# 'fastq_staged' directory with symlinks named per the Space Ranger convention:
#   {sample_name}_S1_L{lane:03d}_{R1|R2}_001.fastq.gz
#
# Input files must end in _<R1/R2/1/2>.fastq.gz.
# Files are grouped by their original sample prefix; each group gets its own lane.
# Prints the sample_name prefix to stdout.

import argparse
import re
import sys
from pathlib import Path

ALLOWED_PATTERN = re.compile(r"^(.+)_(R?[12])\.fastq\.gz$")


def parse_fastq(filename):
    """Return (original_prefix, read_pair) for an allowed FASTQ filename."""
    m = ALLOWED_PATTERN.match(filename)
    if m:
        orig_prefix = m.group(1)
        # Ensure read_pair is in the form R1 or R2
        read_pair = m.group(2) if m.group(2).startswith("R") else f"R{m.group(2)}"
        return orig_prefix, read_pair
    print(
        f"Error: {filename} does not match a recognized allowed FASTQ naming pattern "
        "(_R1.fastq.gz, _R2.fastq.gz, _1.fastq.gz, or _2.fastq.gz).",
        file=sys.stderr,
    )
    sys.exit(1)


def main():
    parser = argparse.ArgumentParser(
        description="Stage and rename allowed FASTQ files for Space Ranger."
    )
    parser.add_argument("fastq_dir", type=Path, help="Directory containing FASTQ files")
    parser.add_argument(
        "sample_name", help="Library ID to use as the shared sample prefix"
    )
    parser.add_argument(
        "staged_dir", type=Path, help="Directory to create with renamed symlinks"
    )
    args = parser.parse_args()

    # Parse args and checks
    fastq_dir = args.fastq_dir
    sample_name = args.sample_name
    staged_dir = args.staged_dir

    if not fastq_dir.is_dir():
        print(
            f"Error: fastq_dir {fastq_dir} does not exist or is not a directory.",
            file=sys.stderr,
        )
        sys.exit(1)
    if not sample_name:
        print("Error: sample_name must not be empty.", file=sys.stderr)
        sys.exit(1)
    if not staged_dir.name:
        print("Error: staged_dir must not be empty.", file=sys.stderr)
        sys.exit(1)

    # List all FASTQ files in order
    # Note their existence was confirmed in nextflow code already
    fastqs = sorted(fastq_dir.glob("*.fastq.gz"))

    # Create directory where we will stage the renamed symlinks
    staged_dir.mkdir(exist_ok=True)

    # Group files by original sample prefix, preserving sort order for lane assignment
    groups = {}
    for f in fastqs:
        orig_prefix, read_pair = parse_fastq(f.name)
        groups.setdefault(orig_prefix, {})[read_pair] = f

    # Create the symlinks with new names according to the Space Ranger convention
    for lane, orig_prefix in enumerate(sorted(groups), start=1):
        for read_pair, f in sorted(groups[orig_prefix].items()):
            new_name = f"{sample_name}_S1_L{lane:03d}_{read_pair}_001.fastq.gz"
            (staged_dir / new_name).symlink_to(f.absolute())

    # send to stdout so we can capture it for input to spaceranger
    print(sample_name)


if __name__ == "__main__":
    main()
