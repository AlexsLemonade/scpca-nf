#!/usr/bin/env python3

# Creates the reference json file for formatting checks for scpca-nf output

import json
from copy import deepcopy

# output reference file
ref_file = "../formatting-reference.json"

# define names of slots and their types for each object type
# --------------------- SCE ---------------------------
# assays
assays = ["counts", "spliced"]
processed_assays = ["logcounts"]
adt_assays = ["counts"]

# cell metadata ------------
cell_metadata = {
    "barcodes": "character",
    "sum": "numeric",
    "detected": "numeric",
    "subsets_mito_sum": "numeric",
    "subsets_mito_detected": "integer",
    "subsets_mito_percent": "numeric",
    "total": "numeric",
}

filtered_cell_metadata = {
    "prob_compromised": "numeric",
    "miQC_pass": "logical",
    "scpca_filter": "character",
    "scDblFinder_class": "character",
    "scDblFinder_score": "numeric",
}

# TODO: don't include these actually because if either of these processes fail then they won't be there
processed_cell_metadata = {
    "sizeFactor": "numeric",
    "cluster": "factor",
}

# row metadata ----------
gene_metadata = {
    "gene_ids": "character",
    "gene_symbol": "character",
    "mean": "numeric",
    "detected": "numeric",
}

# reduced dimensionality ----------
processed_embeddings = ["PCA", "UMAP"]

# alt exps gell/gene metadata ------------

adt_filtered_main_cell_metadata = {"adt_scpca_filter": "character"}

adt_unfiltered_altexp_cell_metadata = {
    "altexps_adt_sum": "numeric",
    "altexps_adt_detected": "integer",
    "altexps_adt_percent": "numeric",
}

adt_unfiltered_altexp_gene_metadata = {
    "adt_id": "character",
    "mean": "numeric",
    "detected": "numeric",
}

# TODO: do we include the colData that's present in filtered projects?
# See https://scpca.readthedocs.io/en/stable/sce_file_contents.html#additional-singlecellexperiment-components-for-cite-seq-libraries-with-adt-tags
# A lot of these are only present in certain conditions, so it doesn't seem worth the check here?

# build base SCE -----------------

unfiltered_sce = {
    "assayNames": assays,
    "colData": {**cell_metadata},
    "rowData": {**gene_metadata},
    "has_submitter": {"colData": {"submitter_celltype_annotation": "character"}},
}

# add in the altExp info
unfiltered_sce["altExp"] = {
    "adt": {
        "assayNames": adt_assays,
        "colData": {**adt_unfiltered_altexp_cell_metadata},
        "rowData": {**adt_unfiltered_altexp_gene_metadata},
    }
}

# build filtered SCE using unfiltered as a base ------
filtered_sce = deepcopy(unfiltered_sce)

filtered_sce["colData"].update({**filtered_cell_metadata})
# use has_adt to specify items that should be in the main object
filtered_sce["has_adt"] = {"colData": {**adt_filtered_main_cell_metadata}}

# build processed SCE with filtered as a base  ------
processed_sce = deepcopy(filtered_sce)
processed_sce["reducedDimNames"] = processed_embeddings

# build and expoort schema --------------
schema = {
    "sce": {
        "unfiltered": unfiltered_sce,
        "filtered": filtered_sce,
        "processed": processed_sce,
    }
    # TODO: Add in anndata
}

with open(ref_file, "w") as f:
    json.dump(schema, f, indent=4)
