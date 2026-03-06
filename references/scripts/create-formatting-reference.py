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
cellhash_assays = ["counts"]

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
    "has_adt": {
        "altexps_adt_sum": "numeric",
        "altexps_adt_detected": "integer",
        "altexps_adt_percent": "numeric",
    },
    "has_cellhash": {
        "altexps_cellhash_sum": "numeric",
        "altexps_cellhash_detected": "integer",
        "altexps_cellhash_percent": "numeric",
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
    # additional columns that are present based on the type of demultiplexing used
    "has_hashedDrops": {"hashedDrops_sampleid": "character"},
    "has_HTODemux": {"HTODemux_sampleid": "character"},
    "has_vireo": {
        "vireo_sampleid": "character",
        "vireo_donor_id": "character",
        "vireo_prob_max": "numeric",
        "vireo_prob_doublet": "numeric",
        "vireo_n_vars": "numeric",
        "vireo_best_singlet": "character",
        "vireo_best_doublet": "character",
        "vireo_doublet_logLikRatio": "numeric",
    },
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

# experiment metadata --------------
unfiltered_experiment_metadata = [
    "library_id",
    "sample_id",
    "project_id",
    "reference_index",
    "total_reads",
    "mapping_tool",
    "tech_version",
    "assay_ontology_term_id",
    "seq_unit",
    "transcript_type",
    "sample_metadata",
    "sample_type",
]

unfiltered_experiment_metadata_conditional = {
    # experiment metadata that's based on mapping tool
    # stored in metadata(sce)$mapping_tool
    "alevin-fry": [
        "salmon_version",
        "mapped_reads",
        "alevinfry_version",
        "af_permit_type",
        "af_resolution",
        "usa_mode",
        "af_num_cells",
        "include_unspliced",
    ],
    "cellranger-multi": [
        "cellranger_version",
        "reference_probeset",
        "pct_mapped_reads",
        "cellranger_num_cells",
    ],
}

filtered_experiment_metadata = unfiltered_experiment_metadata + [
    "filtering_method",
    "miQC_model",
    "prob_compromised_cutoff",
    "scpca_filter_method",
    "min_gene_cutoff",
]

filtered_experiment_metadata_conditional = {
    # if filtering_method is UMI cutoff
    "umi_filtering": ["umi_cutoff"]
}

processed_experiment_metadata = filtered_experiment_metadata + [
    "normalization",
    "highly_variable_genes",
]

processed_experiment_metadata_conditional = {
    **filtered_experiment_metadata_conditional,
    "has_clusters": ["cluster_algorithm", "cluster_weighting", "cluster_nn"],
    # note that celltype_methods exists if either one of these methods is present
    "has_submitter": ["celltype_methods"],
    "has_openscpca": [
        "celltype_methods",
        "openscpca_celltype_module_name",
        "openscpca_celltype_nf_version",
        "openscpca_celltype_release_date",
    ],
    "has_singler": [
        "celltype_methods",
        "singler_results",
        "singler_reference",
        "singler_reference_label",
        "singler_reference_source",
        "singler_reference_version",
    ],
    "has_cellassign": [
        "celltype_methods",
        "cellassign_predictions",
        "cellassign_reference",
        "cellassign_reference_organs",
        "cellassign_reference_source",
        "cellassign_reference_version",
    ],
    "has_scimilarity": [
        "celltype_methods",
        "scimilarity_model",
    ],
    "has_consensus": ["consensus_celltype_methods"],
    "has_infercnv": [
        "infercnv_reference_celltypes",
        "infercnv_num_reference_cells",
        "infercnv_diagnosis_groups",
        "infercnv_status",
        "infercnv_options",
        "infercnv_table",
    ],
    "has_adt": ["adt_scpca_filter_method", "adt_normalization"],
}

# alt exps gell/gene metadata ------------
# adt specific items
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

filtered_altexp_adt_experiment_metadata = ["ambient_profile"]

# cellhash specific items
unfiltered_altexp_cellhash_feature_metadata = {
    "mean": "numeric",
    "detected": "numeric",
}

filtered_altexp_cellhash_feature_metadata = {
    **unfiltered_altexp_cellhash_feature_metadata,
    "barcode_id": "character",
    "sample_id": "character",
}

filtered_altexp_cellhash_cell_metadata_conditional = {
    # indicate columns that are conditionally present based on demultiplexing methods used
    # this is tracked based on demux columns in colData
    # this is the same as what's in the processed object
    "has_hashedDrops": {
        "hashedDrops_Total": "numeric",
        "hashedDrops_Best": "integer",
        "hashedDrops_Second": "integer",
        "hashedDrops_LogFC": "numeric",
        "hashedDrops_LogFC2": "numeric",
        "hashedDrops_Doublet": "logical",
        "hashedDrops_Confident": "logical",
        "hashedDrops_sampleid": "character",
        "hashedDrops_bestsample": "character",
    },
    "has_HTODemux": {
        "HTODemux_maxID": "character",
        "HTODemux_secondID": "character",
        "HTODemux_margin": "numeric",
        "HTODemux_classification": "character",
        "HTODemux_classification.global": "character",
        "HTODemux_hash.ID": "integer",
        "HTODemux_maxsample": "character",
        "HTODemux_sampleid": "character",
    },
}

# build unfiltered SCE -----------------

unfiltered_sce = {
    "assayNames": assays,
    "colData": cell_metadata,
    "rowData": feature_metadata,
    "colData_conditional": cell_metadata_conditional,
    "metadata": unfiltered_experiment_metadata,
    "metadata_conditional": unfiltered_experiment_metadata_conditional,
    "altExp": {
        "adt": {
            "assayNames": adt_assays,
            "rowData": altexp_adt_feature_metadata,
        },
        "cellhash": {
            "assayNames": cellhash_assays,
            "rowData": unfiltered_altexp_cellhash_feature_metadata,
        },
    },
}


# build filtered SCE ------

filtered_sce = {
    "assayNames": assays,
    "colData": filtered_cell_metadata,
    "rowData": feature_metadata,
    "colData_conditional": filtered_cell_metadata_conditional,
    "metadata": filtered_experiment_metadata,
    "metadata_conditional": filtered_experiment_metadata_conditional,
    "altExp": {
        "adt": {
            "assayNames": adt_assays,
            "colData": filtered_altexp_adt_cell_metadata,
            "rowData": altexp_adt_feature_metadata,
            "colData_conditional": filtered_altexp_adt_cell_metadata_conditional,
            "metadata": filtered_altexp_adt_experiment_metadata,
        },
        "cellhash": {
            "assayNames": cellhash_assays,
            "rowData": filtered_altexp_cellhash_feature_metadata,
            "colData_conditional": filtered_altexp_cellhash_cell_metadata_conditional,
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
    "metadata": processed_experiment_metadata,
    "metadata_conditional": processed_experiment_metadata_conditional,
    "altExp": {
        "adt": {
            "assayNames": adt_assays,
            "colData": filtered_altexp_adt_cell_metadata,
            "rowData": altexp_adt_feature_metadata,
            "colData_conditional": processed_altexp_adt_cell_metadata_conditional,
            "metadata": filtered_altexp_adt_experiment_metadata,
        },
        "cellhash": {
            "assayNames": cellhash_assays,
            "rowData": filtered_altexp_cellhash_feature_metadata,
            "colData_conditional": filtered_altexp_cellhash_cell_metadata_conditional,
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
