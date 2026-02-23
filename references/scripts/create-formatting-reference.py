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

cell_metadata_conditional = {
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

filtered_cell_metadata_conditional = {
    **cell_metadata_conditional,
    "has_adt": {"adt_scpca_filter": "character"},
}

# processed cell metadata ------------

processed_cell_metadata_conditional = {
    **filtered_cell_metadata_conditional,
    "has_normalization": {"sizeFactor": "numeric"},
    "has_clusters": {"cluster": "factor"},
    "has_singler": {
        "singler_celltype_annotation": "character",
        "singler_celltype_ontology": "character",
    },
    "has_cellassign": {
        "cellassign_celltype_annotation": "character",
        "cellassign_celltype_ontology": "character",
        "cellassign_max_prediction": "numeric",
    },
    "has_scimilarity": {
        "scimilarity_celltype_annotation": "character",
        "scimilarity_celltype_ontology": "character",
        "scimilarity_min_distance": "numeric",
    },
    "has_consensus": {
        "consensus_celltype_annotation": "character",
        "consensus_celltype_ontology": "character",
    },
    "has_infercnv": {
        "is_infercnv_reference": "logical",
        "infercnv_total_cnv": "integer",
    },
}

# row metadata ----------
feature_metadata = {
    "gene_ids": "character",
    "gene_symbol": "character",
    "mean": "numeric",
    "detected": "numeric",
}

# reduced dimensionality ----------
reduced_dims = ["PCA", "UMAP"]

# alt exps gell/gene metadata ------------
# add main colData columns to conditional metadata
cell_metadata_conditional.update(
    {
        "has_adt": {
            "altexps_adt_sum": "numeric",
            "altexps_adt_detected": "integer",
            "altexps_adt_percent": "numeric",
        },
    }
)

altexp_adt_feature_metadata = {
    "adt_id": "character",
    "mean": "numeric",
    "detected": "numeric",
    "target_type": "character",
}

filtered_altexp_adt_cell_metadata = {
    "zero.ambient": "logical",
    "discard": "logical",
}

filtered_altexp_adt_cell_metadata_conditional = {
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

processed_altexp_adt_cell_metadata_conditional = {
    **filtered_altexp_adt_cell_metadata_conditional,
    "sizeFactor": "numeric",
}

# build unfiltered SCE -----------------

unfiltered_sce = {
    "assayNames": assays,
    "colData": cell_metadata,
    "rowData": feature_metadata,
    "colData_conditional": cell_metadata_conditional,
    "altExp": {
        "adt": {
            "assayNames": adt_assays,
            "rowData": altexp_adt_feature_metadata,
        },
    },
}


# build filtered SCE ------

filtered_sce = {
    "assayNames": assays,
    "colData": filtered_cell_metadata,
    "rowData": feature_metadata,
    "colData_conditional": filtered_cell_metadata_conditional,
    "altExp": {
        "adt": {
            "assayNames": adt_assays,
            "colData": filtered_altexp_adt_cell_metadata,
            "rowData": altexp_adt_feature_metadata,
            "colData_conditional": filtered_altexp_adt_cell_metadata_conditional,
        },
    },
}

# build processed SCE  ------

processed_sce = {
    "assayNames": assays,
    "colData": filtered_cell_metadata,
    "rowData": feature_metadata,
    "colData_conditional": processed_cell_metadata_conditional,
    "reducedDimNames": reduced_dims,
    "altExp": {
        "adt": {
            "assayNames": adt_assays,
            "colData": filtered_altexp_adt_cell_metadata,
            "rowData": altexp_adt_feature_metadata,
            "colData_conditional": processed_altexp_adt_cell_metadata_conditional,
        },
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
