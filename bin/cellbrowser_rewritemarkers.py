#!/usr/bin/env python3

"""
Regenerate per-cluster marker files in a built Cell Browser site.

cbBuild's splitMarkerTable filters out annotated genes because it fails to correctly match
gene symbols to ids, so this script rewrites the marker files using the original annotated marker table and
the features.tsv.gz file to map gene symbols to ids.

Also optionally removes the expression links, as they are not so relevant to these datasets
"""

import argparse
import gzip
import pathlib
import re
import sys

import pandas as pd


def make_cluster_filename(name: str) -> str:
    """
    Sanitize cluster name for use as a filename.
    Mirrors cbBuild's sanitizeName: + to Plus, - to Minus, % to Perc, then strip remaining non-word chars
    """

    # mirrors cbBuild's sanitizeName: + to Plus, - to Minus, % to Perc, then strip remaining non-word chars
    name = name.replace("+", "Plus").replace("-", "Minus").replace("%", "Perc")
    return re.sub(r"[^a-zA-Z_0-9]", "", name)


def make_symbol_dict(features_path: pathlib.Path) -> dict[str, str]:
    """Read features.tsv.gz and return a dict mapping gene symbol to ensembl id."""
    # features.tsv.gz rows:  Ensembl ID, gene symbol
    symbol_dict = {}
    with gzip.open(features_path, "rt") as f:
        for line in f:
            parts = line.rstrip().split("\t")
            if len(parts) >= 2:
                # symbol is the key, id is the value (reverse of features file) for lookup when rewriting markers
                symbol_dict[parts[1]] = parts[0]
    return symbol_dict


def rewrite_markers(
    library_id: str,
    project_id: str,
    cb_data_dir: pathlib.Path,
    outdir: pathlib.Path,
    keep_expression: bool = False,
) -> bool:
    annotated_path = cb_data_dir / project_id / library_id / "markers_annotated.tsv"
    features_path = cb_data_dir / project_id / library_id / "features.tsv.gz"
    markers_out_dir = outdir / project_id / library_id / "markers" / "markers_0"

    if not annotated_path.is_file():
        print(
            f"Skipping {library_id}: no markers_annotated.tsv found at {annotated_path}",
            file=sys.stderr,
        )
        return False
    if not markers_out_dir.is_dir():
        print(
            f"Skipping {library_id}: marker output directory not found at {markers_out_dir}",
            file=sys.stderr,
        )
        return False

    symbol_dict = make_symbol_dict(features_path)

    markers = pd.read_csv(annotated_path, sep="\t")
    cluster_col = markers.columns[0]
    other_cols = [
        c
        for c in markers.columns
        if c not in [cluster_col, "gene", "z_score"]
        # remove expression links (_expr column) unless --keep-expression is set
        and (keep_expression or not c == "_expr")
    ]

    # add id column (Ensembl ID), falling back to the symbol itself if not in the dict
    markers.insert(0, "id", markers["gene"].map(symbol_dict).fillna(markers["gene"]))
    # rename to match cbBuild's expected column names (id, symbol, z_score|float) and preserve other columns
    markers = markers.rename(columns={"gene": "symbol", "z_score": "z_score|float"})
    out_cols = ["id", "symbol", "z_score|float"] + other_cols

    for cluster_name, group in markers.groupby(cluster_col, sort=False):
        out_path = markers_out_dir / (
            make_cluster_filename(str(cluster_name)) + ".tsv.gz"
        )
        group[out_cols].to_csv(out_path, sep="\t", index=False, compression="gzip")

    return True


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--library-ids",
        required=True,
        metavar="ID1,ID2,...",
        help="Comma-separated library IDs to process.",
    )
    parser.add_argument(
        "--project-ids",
        required=True,
        metavar="ID1,ID2,...",
        help="Comma-separated project IDs paired with --library-ids (must be same length).",
    )
    parser.add_argument(
        "--outdir",
        type=pathlib.Path,
        required=True,
        help="Path to the built Cell Browser site directory (cbBuild --outDir).",
    )
    parser.add_argument(
        "--cb-data",
        type=pathlib.Path,
        default=pathlib.Path("cb_data"),
        help="Path to the cb_data staging directory (default: %(default)s).",
    )
    parser.add_argument(
        "--keep-expression",
        action="store_true",
        help="Keep the annotation expression links in the marker files (default removes them).",
    )
    args = parser.parse_args()

    library_ids = args.library_ids.split(",")
    project_ids = args.project_ids.split(",")

    if len(library_ids) != len(project_ids):
        parser.error(
            f"--library-ids ({len(library_ids)}) and --project-ids "
            f"({len(project_ids)}) must have the same number of values."
        )

    success = 0
    for library_id, project_id in zip(library_ids, project_ids):
        if rewrite_markers(
            library_id, project_id, args.cb_data, args.outdir, args.keep_expression
        ):
            success += 1

    print(f"Rewrote marker files for {success}/{len(library_ids)} libraries.")


if __name__ == "__main__":
    main()
