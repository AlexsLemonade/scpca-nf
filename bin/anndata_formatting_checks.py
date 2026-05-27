#!/usr/bin/env python3
# anndata_formatting_checks.py
# Usage: python anndata_formatting_checks.py --anndata_file <file.h5ad>
#            --object_type <unfiltered|filtered|processed>
#            --reference_file <ref.json>
#            [--output_file <output.txt>]
#
# Checks formatting of an AnnData object against a reference JSON.
# Modality (rna or adt) is inferred from the filename.

import argparse
import json
import sys
from pathlib import Path

import anndata
import numpy
import pandas


# Maps reference type strings to Python types for uns value checking.
# The reference uses these strings but actual uns values are numpy scalars/arrays,
# so isinstance is needed rather than comparing type names directly.
_UNS_TYPE_MAP = {
    "str": str,
    "int": (int, numpy.integer),
    "float": (float, numpy.floating),
    "bool": (bool, numpy.bool_),
    "numpy.ndarray": numpy.ndarray,
    "pandas.DataFrame": pandas.DataFrame,
    "NoneType": type(None),
    "dict": dict,
}


def _uns_type_matches(val, expected):
    """Check a uns value against one or more expected type strings."""
    allowed = expected if isinstance(expected, list) else [expected]
    for t in allowed:
        if t is None:
            continue
        # float64 in nested uns (e.g. pca variance/variance_ratio) is a numpy float64 array
        if t == "float64" and isinstance(val, numpy.ndarray):
            if val.dtype == numpy.float64:
                return True
            continue
        types = _UNS_TYPE_MAP.get(t)
        if types is None:
            return True  # unknown type string, skip
        if isinstance(val, types):
            return True
    return False


def check_names_and_types(data, ref, label, slot):
    """
    Check obs or var columns or uns dict for presence and correct type.
    If slot == "uns", data is a dict and types are checked with isinstance.
    If slot == "obs" or "var", data is a DataFrame and types are checked via dtype.name.
    """
    errors = []
    slot_is_uns = slot == "uns"

    for col, expected_type in ref.items():
        if col not in (data.keys() if slot_is_uns else data.columns):
            errors.append(f"Missing '{col}' from {label} {slot}.")

        elif isinstance(expected_type, dict):
            # nested dict: recurse (used for uns keys like "pca")
            errors.extend(check_names_and_types(data[col], expected_type, label, slot))

        elif expected_type is not None:
            if slot_is_uns:
                if not _uns_type_matches(data[col], expected_type):
                    obs_type = type(data[col]).__name__
                    errors.append(
                        f"Type mismatch in '{col}' from {label} {slot}. "
                        f"Expected {expected_type}, but found {obs_type}."
                    )
            else:
                # need to use dtype.name for obs/var columns since they are numpy dtypes
                obs_type = data[col].dtype.name
                if obs_type != expected_type:
                    errors.append(
                        f"Type mismatch in '{col}' from {label} {slot}. "
                        f"Expected {expected_type}, but found {obs_type}."
                    )

    return errors


# for any uns key that is not always present, we want to safely check it's value
# used to set conditionals
def safe_uns_str(uns, key):
    val = uns.get(key)
    return val if isinstance(val, str) else None


# for any uns key (cell type methods) that may be a list or array, we want to safely check and convert to list for conditional checks
def safe_uns_array(uns, key):
    val = uns.get(key)
    if val is None:
        return []
    if isinstance(val, numpy.ndarray):
        return val.tolist()
    if isinstance(val, list):
        return val
    return []


def get_conditionals(adata):
    """
    Build the conditionals dict for an RNA AnnData object.
    Returns a dict where keys are condition names and values are booleans indicating whether the condition is met

    """
    obs_cols = set(adata.obs.columns)
    uns = adata.uns
    # get a list of cell type methods that may be present
    # if not present, then returns [] and all the "has_" conditionals will be False
    celltype_methods = safe_uns_array(uns, "celltype_methods")

    # has adt should be true if any adt-related obs columns are present, or if "altexps_adt_" prefix is used for adt data in uns
    has_adt = (
        any(c.startswith("altexps_adt_") for c in obs_cols)
        or "adt_scpca_filter" in obs_cols
    )
    # is multiplexed if sample_id in uns is a list/array (indicating multiple samples), rather than a single string for non-multiplexed data
    is_multiplexed = isinstance(uns.get("sample_id"), (list, numpy.ndarray))

    # check for submitter annotations based on column and all not being NA
    has_submitter = "submitter" in celltype_methods and (
        "submitter_celltype_annotation" in obs_cols
        and not adata.obs["submitter_celltype_annotation"].isna().all()
    )

    has_neg_ctrl = False
    no_neg_ctrl = False
    if "target_type" in adata.var.columns:
        target_types = adata.var["target_type"]
        has_neg_ctrl = (target_types == "neg_control").any()
        no_neg_ctrl = target_types.isin(["pos_control", "target"]).all()

    return {
        # mapping tools
        "alevin-fry": uns.get("mapping_tool") == "alevin-fry",
        "cellranger-multi": uns.get("mapping_tool") == "cellranger multi",
        # preprocessing steps
        "umi_filtering": uns.get("filtering_method") == "UMI cutoff",
        "has_normalization": "sizeFactor" in obs_cols,
        # clustering and cell type annotation
        "has_clusters": "cluster" in obs_cols,
        "has_consensus": "consensus_celltype_annotation" in obs_cols,
        "has_singler": "singler" in celltype_methods,
        "has_cellassign": "cellassign" in celltype_methods,
        "has_scimilarity": "scimilarity" in celltype_methods,
        "has_openscpca": "openscpca" in celltype_methods,
        "has_submitter": has_submitter,
        "has_infercnv": safe_uns_str(uns, "infercnv_status") == "success",
        # additional modalities
        "has_adt": has_adt,
        "has_cellhash": False,
        "has_hashedDrops": any("hashedDrops" in c for c in obs_cols),
        "has_HTODemux": any("HTODemux" in c for c in obs_cols),
        "has_vireo": any("vireo" in c for c in obs_cols),
        "is_multiplexed": is_multiplexed,
        "has_negative_control": bool(has_neg_ctrl),
        "no_negative_control": bool(no_neg_ctrl),
    }


def check_anndata(adata, ref, label):
    """Run all formatting checks for an AnnData object."""
    errors = []

    # check for contents, raw.X, layers, and obsm keys that are required regardless of conditionals
    if ref.get("has_raw.X") and adata.raw is None:
        errors.append(f"Missing raw.X from {label}.")

    for layer in ref.get("layers", []):
        if layer not in adata.layers:
            errors.append(f"Missing layer '{layer}' from {label}.")

    # get conditionals dict based on adata contents and uns values, to determine which conditional checks to run
    conditions = get_conditionals(adata)
    # only keep true conditions and return as a list of conditionals
    true_conditions = [cond for cond, is_true in conditions.items() if is_true]

    if "obs" in ref:
        errors.extend(check_names_and_types(adata.obs, ref["obs"], label, "obs"))

    if "obs_conditional" in ref:
        for condition, cols_dict in ref["obs_conditional"].items():
            if isinstance(cols_dict, dict) and condition in true_conditions:
                errors.extend(check_names_and_types(adata.obs, cols_dict, label, "obs"))
            elif isinstance(cols_dict, str):
                # condition key is also the column name; always check when present in ref
                errors.extend(
                    check_names_and_types(
                        adata.obs, {condition: cols_dict}, label, "obs"
                    )
                )

    if "var" in ref:
        errors.extend(check_names_and_types(adata.var, ref["var"], label, "var"))

    if "uns" in ref:
        errors.extend(check_names_and_types(adata.uns, ref["uns"], label, "uns"))

    if "uns_conditional" in ref:
        for condition, keys_dict in ref["uns_conditional"].items():
            if isinstance(keys_dict, dict) and condition in true_conditions:
                errors.extend(check_names_and_types(adata.uns, keys_dict, label, "uns"))

    return errors


def main():
    parser = argparse.ArgumentParser(
        description="Check AnnData formatting against a reference JSON."
    )
    parser.add_argument(
        "--anndata_file", required=True, type=Path, help="Path to .h5ad file"
    )
    parser.add_argument(
        "--object_type",
        required=True,
        choices=["unfiltered", "filtered", "processed"],
        help="Type of object: unfiltered, filtered, or processed",
    )
    parser.add_argument(
        "--reference_file", required=True, type=Path, help="Path to reference JSON"
    )
    parser.add_argument(
        "--output_file",
        default="formatting_errors.txt",
        type=Path,
        help=".txt file to save any identified errors",
    )
    args = parser.parse_args()

    if not args.anndata_file.exists():
        sys.exit(f"AnnData file does not exist: {args.anndata_file}")
    if args.output_file.suffix != ".txt":
        sys.exit("output_file must end in .txt")

    if not args.reference_file.exists():
        sys.exit(f"Reference file does not exist: {args.reference_file}")

    # infer modality from filename
    filename = args.anndata_file.stem
    if "_rna" in filename:
        modality = "rna"
    elif "_adt" in filename:
        modality = "adt"
    else:
        sys.exit(
            f"Cannot infer modality (rna or adt) from filename: {args.anndata_file.name}"
        )

    # read in anndata and reference JSON
    adata = anndata.read_h5ad(args.anndata_file)

    with open(args.reference_file) as f:
        all_ref = json.load(f)

    ref = all_ref[args.object_type][modality]
    label = f"{modality.upper()} AnnData"

    # check contents of the anndata object against a reference
    errors = check_anndata(adata, ref, label)

    # export errors to output file, or create empty file if no errors
    output_path = Path(args.output_file)
    if not errors:
        output_path.touch()
    else:
        library_id = adata.uns.get("library_id", "unknown")
        header = (
            f"Formatting errors found for {library_id} {args.object_type} {modality}:"
        )
        with open(output_path, "w") as error_file:
            error_file.write(f"{header}\n\n")
            for error in errors:
                error_file.write(f"{error}\n")
            error_file.write("\n")


if __name__ == "__main__":
    main()
