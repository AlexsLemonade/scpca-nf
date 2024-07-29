#!/usr/bin/env python3
"""
This script takes an AnnData object and checks for the `logcounts` in layers.
If present, `logcounts` is moved to `X` and `X` (which has the raw counts) is moved to `raw.X`

In addition, any DataFrames in `obsm` are conerted to ndarrays, highly variable genes are converted to a `var` column.
If a `pca_meta_file` is provided, PCA variance statistics and standard creation values are in the format expected by scanpy.
"""

import argparse
import os
import re

import anndata
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument(
    "-i",
    "--anndata_file",
    dest="anndata_file",
    required=True,
    help="Path to HDF5 file with processed AnnData object",
)
parser.add_argument(
    "--pca_meta_file",
    dest="pca_meta_file",
    required=False,
    help="Path to file with a table of variance explained by each PCA component",
)
parser.add_argument(
    "--pca_not_centered",
    dest="pca_centered",
    action="store_false",
    help="Indicate that the PCA table is not zero-centered",
)
parser.add_argument(
    "--mask_var",
    dest="mask_var",
    default="highly_variable",
    help="Indicate that variable used for highly variable genes. Use the value 'None' if no highly variable genes were used.",
)
parser.add_argument(
    "-u",
    "--uncompressed",
    dest="compress",
    action="store_false",
    help="Output an uncompressed HDF5 file",
)

args = parser.parse_args()

# compile extension regex
file_ext = re.compile(r"\.(hdf5|h5|h5ad)$", re.IGNORECASE)

# check that input file exists, if it does exist, make sure it's an h5 file
if not os.path.exists(args.anndata_file):
    raise FileExistsError("`input_anndata` does not exist.")
elif not file_ext.search(args.anndata_file):
    raise ValueError(
        "input_anndata must end in either .hdf5, .h5, or .h5ad, and contain a processed AnnData adata."
    )

# if there is a pca_meta_file, check that it exists and has a tsv extension
if args.pca_meta_file:
    if not os.path.exists(args.pca_meta_file):
        raise FileExistsError("`pca_meta_file` does not exist.")
    elif not args.pca_meta_file.casefold().endswith(".tsv"):
        raise ValueError("`pca_meta_file` must end in .tsv.")

# read in anndata
adata = anndata.read_h5ad(args.anndata_file)

# if logcounts is present
if "logcounts" in adata.layers:
    # move X to raw.X by creating the raw adata
    adata.raw = adata
    # move logcounts to X and rename
    adata.X = adata.layers["logcounts"]
    adata.uns["X_name"] = "logcounts"
    del adata.layers["logcounts"]

# convert DataFrames in obsm to ndarrays
for key, value in adata.obsm.items():
    if isinstance(value, pd.DataFrame):
        adata.obsm[key] = value.to_numpy()

# convert highly variable genes to a column
adata.var["highly_variable"] = adata.var.gene_ids.isin(
    adata.uns["highly_variable_genes"]
)

# add pca adata to uns if pca_meta_file is provided in the format created by scanpy
if args.pca_meta_file:
    pca_meta = pd.read_csv(args.pca_meta_file, sep="\t", index_col=0)
    if "variance" not in pca_meta.columns or "variance_ratio" not in pca_meta.columns:
        raise ValueError(
            "`pca_meta_file` must contain columns `variance` and `variance_ratio`"
        )
    pca_object = {
        "param": {
            "zero_center": args.pca_centered,
            "use_highly_variable": args.highly_variable.casefold() != "none",
            "mask_var": args.highly_variable,
        },
        "variance": pca_meta["variance"].to_numpy(),
        "variance_ratio": pca_meta["variance_ratio"].to_numpy(),
    }
    adata.uns["pca"] = pca_object

# export adata
adata.write_h5ad(args.anndata_file, compression="gzip" if args.compress else None)
