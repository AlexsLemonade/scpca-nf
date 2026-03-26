#!/usr/bin/env python3

# Creates the reference json file for formatting checks for scpca-nf output

import json
import copy

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
# used in processed metadata too
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
# this is the same for all object types
feature_metadata = {
    "gene_ids": "character",
    "gene_symbol": "character",
    "mean": "numeric",
    "detected": "numeric",
}

# reduced dimensionality ----------
# only in processed object
reduced_dims = ["PCA", "UMAP"]

# experiment metadata --------------
unfiltered_experiment_metadata = {
    "library_id": "character",
    "sample_id": "character",
    "project_id": "character",
    "reference_index": "character",
    "total_reads": "integer",
    "mapping_tool": "character",
    "tech_version": "character",
    "assay_ontology_term_id": "character",
    "seq_unit": "character",
    "transcript_type": "character",
    "sample_metadata": "data.frame",
    "sample_type": "character",
}

# also used for all adt object types
unfiltered_experiment_metadata_conditional = {
    # experiment metadata that's based on mapping tool
    # stored in metadata(sce)$mapping_tool
    "alevin-fry": {
        "salmon_version": "character",
        "mapped_reads": "integer",
        "alevinfry_version": "character",
        "af_permit_type": "character",
        "af_resolution": "character",
        "usa_mode": "logical",
        "af_num_cells": "integer",
        "include_unspliced": "logical",
    },
    "cellranger-multi": {
        "cellranger_version": "character",
        "reference_probeset": "character",
        "pct_mapped_reads": "numeric",
        "cellranger_num_cells": "integer",
    },
}

filtered_experiment_metadata = {
    **unfiltered_experiment_metadata,
    "filtering_method": "character",
    # this is NA if miQC failed
    "prob_compromised_cutoff": "logical",
    "scpca_filter_method": "character",
    "min_gene_cutoff": "integer",
}

filtered_experiment_metadata_conditional = {
    # if empty drops filtering_method is UMI cutoff
    "umi_filtering": {"umi_cutoff": "numeric"},
    # if miQC is present, these should have specific contents
    # miQC model is an S4 object and is removed from the processed
    "has_miQC": {"miQC_model": "S4", "prob_compromised_cutoff": "numeric"},
}

processed_experiment_metadata = {
    **filtered_experiment_metadata,
    "normalization": "character",
    "highly_variable_genes": "character",
}

processed_experiment_metadata_conditional = {
    ## reuse filtered conditional since miQC model gets removed for processed
    # if empty drops filtering_method is UMI cutoff
    "umi_filtering": {"umi_cutoff": "numeric"},
    # if miQC is present, these should have specific contents
    "has_miQC": {"prob_compromised_cutoff": "numeric"},
    "has_clusters": {
        "cluster_algorithm": "character",
        "cluster_weighting": "character",
        "cluster_nn": "numeric",
    },
    # note that celltype_methods exists if either one of these methods is present
    "has_submitter": {"celltype_methods": "character"},
    "has_openscpca": {
        "celltype_methods": "character",
        "openscpca_celltype_module_name": "character",
        "openscpca_celltype_nf_version": "character",
        "openscpca_celltype_release_date": "character",
    },
    "has_singler": {
        "celltype_methods": "character",
        "singler_results": "DFrame",
        "singler_reference": "character",
        "singler_reference_label": "character",
        "singler_reference_source": "character",
        "singler_reference_version": "character",
        # TODO: these are NA right now because we are using the old SingleR refs
        # only projects that use the new refs will have these set correctly...
        "singler_gene_set_version": "logical",
        "singler_date": "logical",
    },
    "has_cellassign": {
        "celltype_methods": "character",
        "cellassign_predictions": "data.frame",
        "cellassign_reference": "character",
        "cellassign_reference_organs": "character",
        "cellassign_reference_source": "character",
        "cellassign_reference_version": "character",
    },
    "has_scimilarity": {
        "celltype_methods": "character",
        "scimilarity_model": "character",
    },
    "has_consensus": {"consensus_celltype_methods": "character"},
    "has_infercnv": {
        "infercnv_reference_celltypes": "character",
        "infercnv_num_reference_cells": "integer",
        "infercnv_diagnosis_groups": "character",
        "infercnv_status": "character",
        "infercnv_options": "list",
        "infercnv_table": "data.frame",
    },
    "has_adt": {
        "adt_scpca_filter_method": "character",
        "adt_normalization": "character",
    },
}

# alt exps gell/gene metadata ------------
# adt specific items
# same for all object types
altexp_adt_feature_metadata = {
    "adt_id": "character",
    "mean": "numeric",
    "detected": "numeric",
    "target_type": "character",
}

# used for filtered and processed
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

# use for both filtered and processed
filtered_altexp_adt_experiment_metadata = {
    **unfiltered_experiment_metadata,
    "ambient_profile": "character",
}

# cellhash specific items
unfiltered_altexp_cellhash_feature_metadata = {
    "mean": "numeric",
    "detected": "numeric",
}

# used for both filtered and processed
filtered_altexp_cellhash_feature_metadata = {
    **unfiltered_altexp_cellhash_feature_metadata,
    "barcode_id": "character",
    "sample_id": "character",
}

# used for both filtered and processed
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
            "metadata": unfiltered_experiment_metadata,
            "metadata_conditional": unfiltered_experiment_metadata_conditional,
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
            "metadata_conditional": unfiltered_experiment_metadata_conditional,
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
            "metadata_conditional": unfiltered_experiment_metadata_conditional,
        },
        "cellhash": {
            "assayNames": cellhash_assays,
            "rowData": filtered_altexp_cellhash_feature_metadata,
            "colData_conditional": filtered_altexp_cellhash_cell_metadata_conditional,
        },
    },
}

# build and export sce schema --------------
sce_schema = {
    "unfiltered": unfiltered_sce,
    "filtered": filtered_sce,
    "processed": processed_sce,
}

with open(sce_ref_file, "w") as f:
    json.dump(sce_schema, f, indent=4)
    f.write("\n")  # Explicitly add a newline character to appease GitHub

# --------------------- AnnData ---------------------------

# cell and row metadata are dtypes
CELL_ROW_METADATA_MAP = {
    "numeric": "float64",
    "integer": "int32",
    "logical": "bool",
    "character": "category",
}

# outlier types for cell and row metadata
CELL_ROW_METADATA_EXCEPTIONS = {
    "detected": "int32",
    "barcodes": "object",
    "gene_ids": "object",
    "adt_id": "object",
}


def convert_cell_row_metadata_types(metadata):
    for key, value in metadata.items():
        # if the value is a dict, this is conditional metadata and we need to convert the inner values
        # use recursion to do this
        if isinstance(value, dict):
            convert_cell_row_metadata_types(value)
        # check if the key is one of the exceptions where the types aren't what we expect
        elif key in CELL_ROW_METADATA_EXCEPTIONS.keys():
            metadata[key] = CELL_ROW_METADATA_EXCEPTIONS[key]
        # otherwise convert the value as long as the value is in the CELL_ROW_METADATA_MAP
        elif value in CELL_ROW_METADATA_MAP.keys():
            metadata[key] = CELL_ROW_METADATA_MAP[value]

    return metadata


# types for uns metadata entries that are not cell or row metadata
EXPERIMENT_METADATA_MAP = {
    "character": "str",
    "integer": "int",
    "numeric": "float",
    "logical": "bool",
    "data.frame": "pandas.DataFrame",
    "DFrame": "pandas.DataFrame",
}

# outlier types
# most of these stored as numpy.ndarray and not str even though they are considered character in R
EXPERIMENT_METADATA_EXCEPTIONS_MAP = {
    "celltype_methods": "numpy.ndarray",
    "consensus_celltype_methods": "numpy.ndarray",
    "highly_variable_genes": "numpy.ndarray",
    "infercnv_reference_celltypes": "numpy.ndarray",
    "infercnv_options": "numpy.ndarray",
    "transcript_type": "numpy.ndarray",
    "cluster_nn": "int",
}


def convert_experiment_metadata_types(metadata):
    for key, value in metadata.items():
        # account for nested ditionaries in the conditional experiment metadata
        if isinstance(value, dict):
            convert_experiment_metadata_types(value)
        # account for exceptions to the conversion
        elif key in EXPERIMENT_METADATA_EXCEPTIONS_MAP.keys():
            metadata[key] = EXPERIMENT_METADATA_EXCEPTIONS_MAP[key]
        # otherwise convert any values that are present in the map
        elif value in EXPERIMENT_METADATA_MAP.keys():
            metadata[key] = EXPERIMENT_METADATA_MAP[value]

    return metadata


# layer/assays ---------
# used for all objects
layers = ["spliced"]

# unfiltered cell metadata ------------
anndata_specific_obs_metadata = {
    # sample metadata is present in obs for anndata
    # this minimally includes library id and sample id
    "library_id": "category",
    "sample_id": "category",
    # columns that we explicitly add in sce_to_annndata
    "assay_ontology_term_id": "category",
    "suspension_type": "category",
    "is_primary_data": "bool",
}

# if the library is multiplexed, the sample ID doesn't get added
anndata_specific_obs_metadata_conditional = {
    "is_multiplexed": {
        "sample_id": None,
    }
}

obs_metadata = {
    **convert_cell_row_metadata_types(copy.deepcopy(cell_metadata)),
    **anndata_specific_obs_metadata,
}

obs_metadata_conditional = {
    **convert_cell_row_metadata_types(copy.deepcopy(cell_metadata_conditional)),
    **anndata_specific_obs_metadata_conditional,
}

# filtered cell metadata ------------
# used for both filtered and processed
filtered_obs_metadata = {
    **convert_cell_row_metadata_types(copy.deepcopy(filtered_cell_metadata)),
    **anndata_specific_obs_metadata,
}

filtered_cell_metadata_conditional = {
    **convert_cell_row_metadata_types(
        copy.deepcopy(filtered_cell_metadata_conditional)
    ),
    **anndata_specific_obs_metadata_conditional,
}

# processed cell metadata ------------

processed_obs_metadata_conditional = {
    **convert_cell_row_metadata_types(
        copy.deepcopy(processed_cell_metadata_conditional)
    ),
    **anndata_specific_obs_metadata_conditional,
}

# row metadata ----------
# same for all object types
var_metadata = {
    **convert_cell_row_metadata_types(copy.deepcopy(feature_metadata)),
    "feature_is_filtered": "bool",
    "highly_variable": "bool",
}

# reduced dimensionality ----------
processed_obsm = ["X_pca", "X_umap"]

# experiment metadata --------------
anndata_uns_metadata = {
    "X_name": "str",
    "schema_version": "str",
}

unfiltered_uns_metadata = {
    **convert_experiment_metadata_types(copy.deepcopy(unfiltered_experiment_metadata)),
    **anndata_uns_metadata,
}

# also used for all adt object types
unfiltered_uns_metadata_conditional = convert_experiment_metadata_types(
    copy.deepcopy(unfiltered_experiment_metadata_conditional)
)

filtered_uns_metadata = {
    **convert_experiment_metadata_types(copy.deepcopy(filtered_experiment_metadata)),
    **anndata_uns_metadata,
}

filtered_uns_metadata_conditional = {
    # if empty drops filtering_method is UMI cutoff
    "umi_filtering": {"umi_cutoff": "float"},
    # if miQC is present, the prob_compromised_cutoff column should be present
    # no miQC_model since we remove S4 objects
    "has_miQC": {"prob_compromised_cutoff": "float"},
}

processed_uns_metadata = {
    **convert_experiment_metadata_types(copy.deepcopy(processed_experiment_metadata)),
    **anndata_uns_metadata,
    "pca": {
        "param": "dict",
        "variance": "float64",
        "variance_ratio": "float64",
    },
}

processed_uns_metadata_conditional = convert_experiment_metadata_types(
    copy.deepcopy(processed_experiment_metadata_conditional)
)

# alt exps gell/gene metadata ------------
# adt specific items
altexp_adt_var_metadata = convert_cell_row_metadata_types(
    copy.deepcopy(altexp_adt_feature_metadata)
)

# used in both filtered and processed
filtered_altexp_adt_obs_metadata = convert_cell_row_metadata_types(
    copy.deepcopy(filtered_altexp_adt_cell_metadata)
)

filtered_altexp_adt_obs_metadata_conditional = convert_cell_row_metadata_types(
    copy.deepcopy(filtered_altexp_adt_cell_metadata_conditional)
)

processed_altexp_adt_obs_metadata_conditional = convert_cell_row_metadata_types(
    copy.deepcopy(processed_altexp_adt_cell_metadata_conditional)
)

# used in both filtered and processed
filtered_altexp_adt_uns_metadata = convert_experiment_metadata_types(
    copy.deepcopy(filtered_altexp_adt_experiment_metadata)
)

# build unfiltered AnnData -----------------

unfiltered_anndata = {
    "rna": {
        "has_raw.X": False,
        "layers": layers,
        "obs": obs_metadata,
        "var": var_metadata,
        "obs_conditional": obs_metadata_conditional,
        "uns": unfiltered_uns_metadata,
        "uns_conditional": unfiltered_uns_metadata_conditional,
    },
    "adt": {
        "var": altexp_adt_var_metadata,
    },
}


# build filtered anndata ------

filtered_anndata = {
    "rna": {
        "has_raw.X": False,
        "layers": layers,
        "obs": filtered_obs_metadata,
        "var": var_metadata,
        "obs_conditional": filtered_cell_metadata_conditional,
        "uns": filtered_uns_metadata,
        "uns_conditional": filtered_uns_metadata_conditional,
    },
    "adt": {
        "obs": filtered_altexp_adt_obs_metadata,
        "var": altexp_adt_var_metadata,
        "obs_conditional": filtered_altexp_adt_obs_metadata_conditional,
        "uns": filtered_altexp_adt_uns_metadata,
    },
}

# build processed AnnData  ------

processed_anndata = {
    "rna": {
        "has_raw.X": True,
        "layers": layers,
        "obs": filtered_obs_metadata,
        "var": var_metadata,
        "obs_conditional": processed_obs_metadata_conditional,
        "obsm": processed_obsm,
        "uns": processed_uns_metadata,
        "uns_conditional": processed_uns_metadata_conditional,
    },
    "adt": {
        "obs": filtered_altexp_adt_obs_metadata,
        "var": altexp_adt_var_metadata,
        "obs_conditional": processed_altexp_adt_obs_metadata_conditional,
        "uns": filtered_altexp_adt_uns_metadata,
    },
}

# build and export anndata schema --------------
anndata_schema = {
    "unfiltered": unfiltered_anndata,
    "filtered": filtered_anndata,
    "processed": processed_anndata,
}

with open(anndata_ref_file, "w") as f:
    json.dump(anndata_schema, f, indent=4)
    f.write("\n")  # Explicitly add a newline character to appease GitHub
