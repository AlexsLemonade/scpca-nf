#!/usr/bin/env python3

# Script to annotate processed ScPCA objects with SCimilarity
# Follows this tutorial: https://genentech.github.io/scimilarity/notebooks/cell_annotation_tutorial.html

# Adapted from: https://github.com/AlexsLemonade/OpenScPCA-nf/blob/v0.1.6/modules/cell-type-scimilarity/resources/usr/bin/run-scimilarity.py

import argparse
import sys
from pathlib import Path
import random

import anndata
import pandas
from scipy.sparse import csr_matrix
from scimilarity import CellAnnotation
from scimilarity.utils import align_dataset, lognorm_counts


def format_scimilarity(adata: anndata.AnnData) -> anndata.AnnData:
    """
    Creates a new AnnData object formatted for running SCimilarity:

    - Gene symbols, taken from the 'gene_symbol' column in 'adata.var' are used as var_names
    If any duplicate gene symbols are found, they are collapsed by summing counts.
    - The summed counts matrix is then stored in the 'counts' layer of the new AnnData object.

    Parameters
    ----------
    adata : AnnData
        Input AnnData object, with "gene_symbol" column in `adata.var`.

    Returns
    -------
    AnnData
        New AnnData object with gene symbols as var_names and counts stored in the 'counts' layer.
    """

    # Check that gene symbol column exists
    if "gene_symbol" not in adata.var.columns:
        raise ValueError(
            "The input AnnData object must have a 'gene_symbol' column in `adata.var`."
        )

    # set gene symbols as var_names and make sure X has the raw counts
    adata.var_names = adata.var["gene_symbol"].astype(str)
    adata.X = adata.raw.X

    # create a DataFrame with raw counts, dropping anything that doesn't have a gene symbol
    counts_df = adata.to_df().drop(columns=["nan"])
    # Collapse duplicates by summing and make sparse
    collapsed_df = counts_df.T.groupby(level=0).sum().T

    # Build new AnnData with collapsed counts stored as layers
    # this is expected by SCimilarity
    adata_collapsed = anndata.AnnData(
        obs=pandas.DataFrame(index=collapsed_df.index),  # original cell barcodes
        var=pandas.DataFrame(
            index=collapsed_df.columns
        ),  # gene symbols after collapsing
        layers={"counts": csr_matrix(collapsed_df)},
    )

    return adata_collapsed


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Annotate cells using SCimilarity",
    )
    parser.add_argument(
        "--model_dir",
        type=Path,
        required=True,
        help="Path to the directory containing the SCimilarity foundational model",
    )
    parser.add_argument(
        "--processed_h5ad_file",
        type=Path,
        required=True,
        help="Path to the processed AnnData object stored as an h5ad file",
    )
    parser.add_argument(
        "--ontology_map_file",
        type=Path,
        required=True,
        help="Path to TSV file containing Cell Ontology identifiers for annotation terms",
    )
    parser.add_argument(
        "--predictions_tsv",
        type=Path,
        required=True,
        help="Path to the output TSV file with cell type annotations",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=2025,
        help="Random seed to ensure reproducibility",
    )
    arg = parser.parse_args()

    # Set the seed
    random.seed(arg.seed)

    ################################################
    ########### Input argument checks ##############
    ################################################
    arg_error = False

    # Check that input files exist
    if not arg.model_dir.is_dir():
        print(
            f"The provided reference SCimilarity model could not be found at: {arg.model_dir}.",
            file=sys.stderr,
        )
        arg_error = True
    if not arg.processed_h5ad_file.is_file():
        print(
            f"The provided input H5AD file could not be found at: {arg.processed_h5ad_file}.",
            file=sys.stderr,
        )
        arg_error = True
    if not arg.ontology_map_file.is_file():
        print(
            f"The ontology map file could not be found at: {arg.ontology_map_file}.",
            file=sys.stderr,
        )
        arg_error = True

    # Exit if error(s)
    if arg_error:
        sys.exit(1)

    ################################################
    ################ Prep input data ###############
    ################################################

    # Read in model
    scimilarity_model = CellAnnotation(model_path=arg.model_dir)

    # read in ontology identifiers
    ontology_map = pandas.read_csv(
        arg.ontology_map_file, sep="\t", index_col="scimilarity_celltype_annotation"
    )

    # Read and make sure object formatting is correct
    processed_anndata = anndata.read_h5ad(arg.processed_h5ad_file)
    processed_anndata = format_scimilarity(processed_anndata)

    # Preprocess the data
    # Align the query dataset to the reference model
    processed_anndata = align_dataset(processed_anndata, scimilarity_model.gene_order)
    # Log-normalize the counts
    processed_anndata = lognorm_counts(processed_anndata)

    ################################################
    ############### Run Scimilarity ###############
    ################################################

    # compute embeddings
    processed_anndata.obsm["X_scimilarity"] = scimilarity_model.get_embeddings(
        processed_anndata.X
    )

    # Predict cell types
    predictions, nn_idxs, nn_dists, nn_stats = scimilarity_model.get_predictions_knn(
        processed_anndata.obsm["X_scimilarity"], weighting=True
    )

    ################################################
    ################ Export annotations ############
    ################################################

    # prepare the predictions with min distance for export
    predictions_df = pandas.DataFrame(
        {
            "barcode": processed_anndata.obs_names.to_list(),
            "scimilarity_celltype_annotation": predictions.values,
            "min_dist": nn_stats["min_dist"],
        }
    )
    # add in ontology IDs
    predictions_df = predictions_df.join(
        ontology_map, on="scimilarity_celltype_annotation", validate="many_to_one"
    )

    # export TSV
    predictions_df.to_csv(arg.predictions_tsv, sep="\t", index=False)


if __name__ == "__main__":
    main()
