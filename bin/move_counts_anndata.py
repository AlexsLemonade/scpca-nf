#!/usr/bin/env python3

# This script takes an AnnData object and checks for the `logcounts`
# in layers. If present, `logcounts` is moved to `X` and `X` (which has the raw counts)
# is moved to `raw.X`

import os
import anndata as adata
import argparse
import re

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input_hdf5_file',
                    dest = 'input_hdf5_file',
                    required = True,
                    help = 'Path to HDF5 file(s) with processed AnnData object. If more than one file, separate with comma.')

args = parser.parse_args()

# separate input files
input_files = args.input_hdf5_file.split(",")

# compile extension regex
file_ext = re.compile(r"\.hdf5$|.h5$", re.IGNORECASE)

for file in input_files:
    # check that input file exists, if it does exist, make sure it's an h5 file
    if not os.path.exists(file):
        raise FileExistsError("Make sure that all input_hdf5_file exist.")
    elif not file_ext.search(file):
        raise ValueError("--input_hdf5_file must end in either .hdf5 or .h5 and contain a processed AnnData object.")

    # read in anndata
    object = adata.read_h5ad(file)

    # if logcounts is present
    if "logcounts" in object.layers.keys():
        # move X to raw.X by creating the raw object
        object.raw = object

        # move logcounts to X and rename
        object.X = object.layers["logcounts"]
        object.uns["X_name"] = "logcounts"

        # export object
        object.write_h5ad(file)
