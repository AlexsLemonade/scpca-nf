#!/usr/bin/env python3

# This script takes an AnnData object and a reference marker gene matrix
# where genes are rows and columns are cell types
# CellAssign is run and the predictions matrix is returned


import argparse
import os
import re

import anndata as adata
import numpy as np
import pandas as pd
import scvi
from scvi.external import CellAssign

parser = argparse.ArgumentParser()
parser.add_argument(
    "-i",
    "--anndata_file",
    dest="anndata_file",
    required=True,
    help="Path to HDF5 file with processed AnnData object to annotate",
)
parser.add_argument(
    "-o",
    "--output_predictions",
    dest="output_predictions",
    required=True,
    help="Path to file to save the predictions, should end in tsv",
)
parser.add_argument(
    "-r",
    "--reference",
    dest="reference",
    required=True,
    help="Path to marker by cell type reference file, should end in tsv",
)
parser.add_argument(
    "-s",
    "--seed",
    dest="seed",
    type=int,
    default=None,
    help="Random seed to set for scvi.",
)
parser.add_argument(
    "-t",
    "--threads",
    dest="threads",
    default=1,
    type=int,
    help="Number of threads to use when training the CellAssign model",
)

args = parser.parse_args()

# set seed
scvi.settings.seed = args.seed

# define the number of threads to use
scvi.settings.num_threads = args.threads

# compile extension regex
file_ext = re.compile(r"\.(hdf5|h5|h5ad)$", re.IGNORECASE)

# check that input file exists, if it does exist, make sure it's an h5 file
if not os.path.exists(args.anndata_file):
    raise FileExistsError("--anndata_file file not found.")
elif not file_ext.search(args.anndata_file):
    raise ValueError(
        "--anndata_file must end in .h5ad, .hdf5, or .h5 and contain a processed AnnData object."
    )

# check that marker file exists and make sure its a tsv
if not os.path.exists(args.reference):
    raise FileExistsError("--reference file not found.")

# make sure output file path is tsv file
if not args.output_predictions.casefold().endswith(".tsv"):
    raise ValueError("--output_predictions must provide a file path ending in tsv")

# read in references as marker gene tables
ref_matrix = pd.read_csv(args.reference, sep="\t", index_col="ensembl_id")

# file path to annotated sce
annotated_adata = adata.read_h5ad(args.anndata_file)

# subset anndata to contain only genes in the reference file
# note that the gene names must be the rownames of the reference matrix
# first get a list of shared genes
shared_genes = list(set(ref_matrix.index) & set(annotated_adata.var_names))

# check that shared_genes actually has some genes
if not shared_genes:
    raise ValueError(
        "--reference does not include any genes found in the provided --anndata_file."
    )

# create a new anndata object with only shared genes
subset_adata = annotated_adata[:, shared_genes].copy()
subset_adata.X = subset_adata.X.tocsr()

# add size factor to subset adata (calculated from full data)
lib_size = annotated_adata.X.sum(1)
subset_adata.obs["size_factor"] = lib_size / np.mean(lib_size)

# only run CellAssign if enough cells
if subset_adata.n_obs < 30:
    # make a predictions file that just has the barcode column
    barcodes_column = subset_adata.obs_names.to_list()
    predictions = pd.DataFrame(barcodes_column, columns=["barcode"])
else:
    # train and assign cell types
    scvi.external.CellAssign.setup_anndata(subset_adata, size_factor_key="size_factor")
    model = CellAssign(subset_adata, ref_matrix)
    model.train()
    predictions = model.predict()
    predictions["barcode"] = subset_adata.obs_names


# write out predictions as tsv
predictions.to_csv(args.output_predictions, sep="\t", index=False)
