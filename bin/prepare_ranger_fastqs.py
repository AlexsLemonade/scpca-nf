#!/usr/bin/env python3
# prepare_ranger_fastqs.py
# Usage: prepare_ranger_fastqs.py --fastq-dir <fastq_dir> --staged-dir <staged_dir>
#
# Stages FASTQ files for Space Ranger or Cell Ranger input by creating a directory
# of symlinks and captures the sample prefix(es) for input.
#
# Conformant files (already 10x-formatted) are symlinked with their original names.
#
# Non-conformant allowed files (_R1/R2.fastq.gz or _1/2.fastq.gz, optionally with
# a _\d{3} suffix) are symlinked with names following the 10x Genomics convention:
#   {prefix}_S1_{R1|R2}_001.fastq.gz
#
# The comma-separated list of sample prefixes is printed to stdout.
# Exits non-zero if any file does not match either recognized pattern.

import argparse
import re
import sys
from pathlib import Path

CONFORMANT_PATTERN = re.compile(r"^(.+)_S\d+_(L\d+_)?[RI][12]_\d{3}\.fastq\.gz$")
ALLOWED_PATTERN = re.compile(r"^(?P<prefix>.+)_(?P<read>[RI]?[12])(?:_\d{3})?\.fastq\.gz$")


def parse_allowed_fastq(filename):
    """Return (original_prefix, read_pair) for an allowed non-conformant FASTQ filename."""
    m = ALLOWED_PATTERN.match(filename)
    if m:
        orig_prefix = m.group("prefix")
        read = m.group("read")
        read_pair = read if read.startswith("R") else f"R{read}"
        return orig_prefix, read_pair
    print(
        f"Error: {filename} does not match a recognized allowed FASTQ naming pattern "
        "(_R1.fastq.gz, _R2.fastq.gz, _1.fastq.gz, or _2.fastq.gz).",
        file=sys.stderr,
    )
    sys.exit(1)


def main():
    parser = argparse.ArgumentParser(
        description="Stage FASTQ files for input to Space Ranger or Cell Ranger."
    )
    parser.add_argument(
        "--fastq-dir", required=True, type=Path, help="Directory containing FASTQ files"
    )
    parser.add_argument(
        "--staged-dir",
        required=True,
        type=Path,
        help="Directory to create with symlinks",
    )
    args = parser.parse_args()

    if not args.fastq_dir.is_dir():
        print(
            f"Error: fastq-dir {args.fastq_dir} does not exist or is not a directory.",
            file=sys.stderr,
        )
        sys.exit(1)

    # List all FASTQ files
    # Note their existence was confirmed in nextflow code already
    fastqs = list(args.fastq_dir.glob("*.fastq.gz"))

    # Create directory where we will stage all symlinks
    args.staged_dir.mkdir()

    return_prefixes = set()
    for f in fastqs:
        # don't change the name if it's already conformant
        if m := CONFORMANT_PATTERN.match(f.name):
            prefix = m.group(1)
            new_name = f.name
        # add _S1_ and _001 to conform to 10x convention
        else:
            prefix, read_pair = parse_allowed_fastq(f.name)
            new_name = f"{prefix}_S1_{read_pair}_001.fastq.gz"
        (args.staged_dir / new_name).symlink_to(f.absolute())
        return_prefixes.add(prefix)

    if not return_prefixes:
        print(
            "Error: Could not identify sample prefixes from the FASTQ files.",
            file=sys.stderr,
        )
        sys.exit(1)

    # Send to stdout so we can capture all prefixes for input to Space Ranger/Cell Ranger
    print(",".join(return_prefixes))


if __name__ == "__main__":
    main()
