#!/usr/bin/env python3

# This script takes an AnnData object and a reference marker gene matrix
# where genes are rows and columns are cell types
# CellAssign is run and the predictions matrix is returned


import os
import anndata as adata
import scvi
from scvi.external import CellAssign
from scvi import settings
import pandas as pd
import numpy as np
import argparse
import re

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input_anndata',
                    dest = 'input_anndata',
                    required = True,
                    help = 'Path to HDF5 file with processed AnnData object to annotate')
parser.add_argument('-o', '--output_predictions',
                    dest = 'output_predictions',
                    required = True,
                    help = 'Path to directory to save the predictions, should end in tsv')
parser.add_argument('-r', '--reference',
                    dest = 'reference',
                    required = True,
                    help = 'Path to marker by cell type reference file')
parser.add_argument('-s', '--seed',
                    dest = 'seed',
                    type=int,
                    default = None,
                    help = 'Random seed to set for scvi.')
parser.add_argument('-t', '--threads',
                    dest = 'threads',
                    default = 1,
                    type = int,
                    help = 'Number of threads to use when training the CellAssign model')

args = parser.parse_args()

# set seed
scvi.settings.seed = args.seed

# define the number of threads to use
settings.num_threads = args.threads

# compile extension regex
file_ext = re.compile(r"\.hdf5$|.h5$", re.IGNORECASE)

# check that input file exists, if it does exist, make sure it's an h5 file
if not os.path.exists(args.input_anndata):
    raise FileExistsError("--input_anndata file not found.")
elif not file_ext.search(args.input_anndata):
    raise ValueError("--input_anndata must end in either .hdf5 or .h5 and contain a processed AnnData object.")

# check that marker file exists and make sure its a csv
if not os.path.exists(args.reference):
    raise FileExistsError("--reference file not found.")
elif not ".csv" in args.reference:
    raise ValueError("--reference must be a csv file")

# make sure output file path is tsv file
if not ".tsv" in args.output_predictions:
    raise ValueError("--output_predictions must provide a file path ending in tsv")

# check that output file directory exists and create directory if doesn't exist
predictions_dir = os.path.dirname(args.output_predictions)
os.makedirs(predictions_dir, exist_ok = True)

# read in references as marker gene tables
ref_matrix = pd.read_csv(args.reference, sep = "\t", index_col='ensembl_id')

# file path to annotated sce
annotated_adata = adata.read_h5ad(args.input_anndata)

# subset anndata to contain only genes in the reference file
# note that the gene names must be the rownames of the reference matrix
subset_adata = annotated_adata[:, ref_matrix.index].copy()
subset_adata.X = subset_adata.X.tocsr()

# add size factor to subset adata (calculated from full data)
lib_size = annotated_adata.X.sum(1)
subset_adata.obs["size_factor"] = lib_size / np.mean(lib_size)

# train and assign cell types
scvi.external.CellAssign.setup_anndata(subset_adata, size_factor_key="size_factor")
model = CellAssign(subset_adata, ref_matrix)
model.train()
predictions = model.predict()
predictions['barcode'] = subset_adata.obs_names

# write out predictions as tsv
predictions.to_csv(args.output_predictions, sep = "\t", index = False)
