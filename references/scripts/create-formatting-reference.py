#!/usr/bin/env python3

# Creates the reference json file for formatting checks for scpca-nf output

import json

# output reference files
sce_ref_file = "../sce-formatting-reference.json"
anndata_ref_file = "../anndata-formatting-reference.json"

# define names of slots and their types for each object type
# --------------------- SCE ---------------------------
# assays
assays = ["counts", "spliced"]
processed_assays = ["logcounts"]
adt_assays = ["counts"]

# unfiltered cell metadata ------------
cell_metadata = {
    "barcodes": "character",
    "sum": "numeric",
    "detected": "numeric",
    "subsets_mito_sum": "numeric",
    "subsets_mito_detected": "integer",
    "subsets_mito_percent": "numeric",
    "total": "numeric",
}

conditional_cell_metadata = {
    "has_submitter": {"submitter_celltype_annotation": "character"},
    "has_openscpca": {
        "openscpca_celltype_annotation": "character",
        "openscpca_celltype_ontology": "character",
    },
}

# filtered cell metadata ------------
filtered_cell_metadata = {
    **cell_metadata,
    "prob_compromised": "numeric",
    "miQC_pass": "logical",
    "scpca_filter": "character",
    "scDblFinder_class": "character",
    "scDblFinder_score": "numeric",
}

filtered_conditional_cell_metadata = {
    **conditional_cell_metadata,
    "has_adt": {"adt_scpca_filter": "character"},
}

# processed cell metadata ------------

processed_conditional_cell_metadata = {
    **filtered_conditional_cell_metadata,
    "has_normalization": {"sizeFactor": "numeric"},
    "has_clusters": {"cluster": "factor"},
    # TODO: should add another level of nesting with has_singler, has_cellassign, etc?
    "has_celltyping": {
        "singler_celltype_annotation": "character",
        "singler_celltype_ontology": "character",
        "cellassign_celltype_annotation": "character",
        "cellassign_celltype_ontology": "character",
        "cellassign_max_prediction": "numeric",
        "scimilarity_celltype_annotation": "character",
        "scimilarity_celltype_ontology": "character",
        "scimilarity_min_distance": "numeric",
        "consensus_celltype_annotation": "character",
        "consensus_celltype_ontology": "character",
    },
    "has_infercnv": {
        "is_infercnv_reference": "logical",
        "infercnv_total_cnv": "integer",
    },
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

adt_unfiltered_altexp_cell_metadata = {
    "altexps_adt_sum": "numeric",
    "altexps_adt_detected": "integer",
    "altexps_adt_percent": "numeric",
}

adt_unfiltered_altexp_gene_metadata = {
    "adt_id": "character",
    "mean": "numeric",
    "detected": "numeric",
    "target_type": "character",
}

adt_filtered_altexp_cell_metadata = {
    **adt_unfiltered_altexp_cell_metadata,
    "zero.ambient": "logical",
    "discard": "logical",
}

adt_filtered_conditional_altexp_cell_metadata = {
    # indicate columns that are conditionally present based on negative control
    # this is tracked by target column of rowData
    "has_negative_control": {
        "sum.controls": "integer",
        "high.control": "logical",
    },
    "no_negative_control": {
        "ambient.scale": "numeric",
        "high.ambient": "logical",
    },
}

adt_processed_conditional_altexp_cell_metadata = {
    **adt_filtered_conditional_altexp_cell_metadata,
    "sizeFactor": "numeric",
}

# build unfiltered SCE -----------------

unfiltered_sce = {
    "assayNames": assays,
    "colData": cell_metadata,
    "rowData": gene_metadata,
    "colData_conditional": conditional_cell_metadata,
    "altExp": {
        "adt": {
            "assayNames": adt_assays,
            "colData": adt_unfiltered_altexp_cell_metadata,
            "rowData": adt_unfiltered_altexp_gene_metadata,
        }
    },
}


# build filtered SCE ------

filtered_sce = {
    "assayNames": assays,
    "colData": filtered_cell_metadata,
    "rowData": gene_metadata,
    "colData_conditional": filtered_conditional_cell_metadata,
    "altExp": {
        "adt": {
            "assayNames": adt_assays,
            "colData": adt_filtered_altexp_cell_metadata,
            "rowData": adt_unfiltered_altexp_gene_metadata,
            "colData_conditional": adt_filtered_conditional_altexp_cell_metadata,
        }
    },
}

# build processed SCE  ------

processed_sce = {
    "assayNames": assays,
    "colData": filtered_cell_metadata,
    "rowData": gene_metadata,
    "colData_conditional": processed_conditional_cell_metadata,
    "reducedDimNames": processed_embeddings,
    "altExp": {
        "adt": {
            "assayNames": adt_assays,
            "colData": adt_filtered_altexp_cell_metadata,
            "rowData": adt_unfiltered_altexp_gene_metadata,
            "colData_conditional": adt_processed_conditional_altexp_cell_metadata,
        }
    },
}

# build and expoort sce schema --------------
sce_schema = {
    "unfiltered": unfiltered_sce,
    "filtered": filtered_sce,
    "processed": processed_sce,
}

with open(sce_ref_file, "w") as f:
    json.dump(sce_schema, f, indent=4)
    f.write("\n")  # Explicitly add a newline character to appease GitHub
