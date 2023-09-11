#!/usr/bin/env python3

# This script takes an AnnData object and checks for the `logcounts`
# in layers. If present, `logcounts` is moved to `X` and `X` (which has the raw counts)
# is moved to `raw.X`

import argparse
import os
import re

import anndata as adata

parser = argparse.ArgumentParser()
parser.add_argument(
    "-i",
    "--anndata_file",
    dest="anndata_file",
    required=True,
    help="Path to HDF5 file with processed AnnData object",
)

args = parser.parse_args()

# compile extension regex
file_ext = re.compile(r"\.hdf5$|.h5$", re.IGNORECASE)

# check that input file exists, if it does exist, make sure it's an h5 file
if not os.path.exists(args.anndata_file):
    raise FileExistsError("`input_anndata` does not exist.")
elif not file_ext.search(args.anndata_file):
    raise ValueError(
        "--input_anndata must end in either .hdf5 or .h5 and contain a processed AnnData object."
    )

# read in anndata
object = adata.read_h5ad(args.anndata_file)

# if logcounts is present
if "logcounts" in object.layers:
    # move X to raw.X by creating the raw object
    object.raw = object
    # move logcounts to X and rename
    object.X = object.layers["logcounts"]
    object.uns["X_name"] = "logcounts"

    # add schema to uns
    object.uns["schema_version"] = "3.0.0"

    # export object
    object.write_h5ad(args.anndata_file)
