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
import numpy as np
import pandas as pd


# Maps simplified type names to pandas dtype-checking functions.
# These work for both numpy arrays and pandas Series regardless of backend.
_DTYPE_CHECKERS = {
    "int": pd.api.types.is_integer_dtype,
    "float": pd.api.types.is_float_dtype,
    "string": pd.api.types.is_string_dtype,
    "bool": pd.api.types.is_bool_dtype,
    # account for deprecated pd.CategoricalDtype in older pandas versions
    "category": lambda obj: isinstance(getattr(obj, "dtype", obj), pd.CategoricalDtype),
}

_BUILTIN_TYPES = {
    "str": str,
    "int": (int, np.integer),  # covers int64, int32, etc.
    "float": (float, np.floating),  # covers float64, float32, etc.
    "bool": (bool, np.bool_),  # covers numpy.bool_
    "NoneType": type(None),
    "dict": dict,
    "numpy.ndarray": np.ndarray,
}


def _check_type(obj, expected_type):
    """
    Check whether obj matches expected_type.

    For pandas Series and numpy arrays, uses pd.api.types functions so that
    type checking is robust across pandas backends and numpy dtype variations.
    For plain Python objects, falls back to isinstance with built-in types.

    Returns True if the type matches, False otherwise.
    """

    # first check for numpy arrays, since these can be present in uns
    # and require a different check than plain Python objects or pandas Series
    if expected_type == "numpy.ndarray":
        return isinstance(obj, np.ndarray)

    # if the object is a Series or array, use the dtype checker
    # this runs pd.api.types to check the dtype
    if isinstance(obj, (pd.Series, np.ndarray)):
        checker = _DTYPE_CHECKERS.get(expected_type)
        if checker is None:
            return False
        return checker(obj)
    else:
        # otherwise check against built in python or numpy types
        # account for expected_type being a single type or a comma-separated list of types in the reference
        allowed = (
            expected_type
            if isinstance(expected_type, list)
            else [t.strip() for t in str(expected_type).split(",")]
        )
        py_types = [_BUILTIN_TYPES.get(t) for t in allowed]
        return any(isinstance(obj, t) for t in py_types if t is not None)


def check_names_and_types(data, ref, label, slot):
    """
    Check obs/var columns or uns dict for presence and correct type.

    For obs/var slots, data is a DataFrame; column values are pd.Series.
    For uns slot, data is a dict; values may be plain Python objects,
    numpy arrays, or pandas Series — _check_type handles all cases.
    """
    errors = []
    keys = data.keys() if hasattr(data, "keys") else data.columns

    for key, expected_type in ref.items():
        if key not in keys:
            errors.append(f"Missing '{key}' from {label} {slot}.")

        elif isinstance(expected_type, dict):
            # nested dict: recurse (used for uns keys like "pca")
            errors.extend(
                check_names_and_types(data[key], expected_type, label, f"{slot}:{key}")
            )

        elif expected_type is not None:
            obj = data[key]
            if not _check_type(obj, expected_type):
                # get the observed type based on if the object is a Series/array or a plain Python object
                obs_type = (
                    obj.dtype.name
                    if isinstance(obj, (pd.Series, np.ndarray))
                    else type(obj).__name__
                )
                errors.append(
                    f"Type mismatch in '{key}' from {label} {slot}. "
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
    if isinstance(val, np.ndarray):
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
    is_multiplexed = isinstance(uns.get("sample_id"), (list, np.ndarray))

    # check for submitter annotations based on column and all not being NA
    has_submitter = "submitter" in celltype_methods and (
        "submitter_celltype_annotation" in obs_cols
        and not adata.obs["submitter_celltype_annotation"].isna().all()
    )

    # set neg post control based on target types if adt
    has_neg_ctrl = False
    no_neg_ctrl = False
    if "target_type" in adata.var.columns:
        target_types = adata.var["target_type"]
        has_neg_ctrl = (target_types == "neg_control").any()
        no_neg_ctrl = target_types.isin(["pos_control", "target"]).all()

    condition_tests = {
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

    return set(cond for cond, is_true in condition_tests.items() if is_true)


def expand_merged_var_ref(var, ref_var):
    """
    Expand the var reference for merged AnnData objects to include per-library columns.

    For merged objects, var contains `libraryid.mean` and `libraryid.detected`
    columns for each library.
    This function finds those columns and adds them to the reference
    with the same expected types as `mean` and `detected`.
    """
    expanded = ref_var.copy()
    
    # get mean and detected values, while removing them so they don't get flagged as missing
    mean_type = expanded.pop("mean")
    detected_type = expanded.pop("detected")

    for col in var.columns:
        if col.endswith(".mean"):
            expanded[col] = mean_type
        elif col.endswith(".detected"):
            expanded[col] = detected_type

    return expanded


def check_anndata(adata, ref, label, object_type):
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

    if "obs" in ref:
        errors.extend(check_names_and_types(adata.obs, ref["obs"], label, "obs"))

    if "obs_conditional" in ref:
        for condition, cols_dict in ref["obs_conditional"].items():
            if isinstance(cols_dict, dict) and condition in conditions:
                errors.extend(check_names_and_types(adata.obs, cols_dict, label, "obs"))
            elif isinstance(cols_dict, str):
                # condition key is also the column name; always check when present in ref
                errors.extend(
                    check_names_and_types(
                        adata.obs, {condition: cols_dict}, label, "obs"
                    )
                )

    if "var" in ref:
        if object_type == "merged":
            ref["var"] = expand_merged_var_ref(adata.var, ref["var"])
        errors.extend(check_names_and_types(adata.var, ref["var"], label, "var"))

    if "uns" in ref:
        errors.extend(check_names_and_types(adata.uns, ref["uns"], label, "uns"))

    if "uns_conditional" in ref:
        for condition, keys_dict in ref["uns_conditional"].items():
            if isinstance(keys_dict, dict) and condition in conditions:
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
        "--object_id",
        required=True,
        help="Identifier for the object (e.g. library ID) to include in error messages",
    )
    parser.add_argument(
        "--object_type",
        required=True,
        choices=["unfiltered", "filtered", "processed", "merged"],
        help="Type of object: unfiltered, filtered, processed, or merged",
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
    errors = check_anndata(adata, ref, label, args.object_type)

    # export errors to output file, or create empty file if no errors
    output_path = Path(args.output_file)
    if not errors:
        output_path.touch()
    else:
        header = f"Formatting errors found for {args.object_id} {args.object_type} {modality} AnnData object:"
        with open(output_path, "w") as error_file:
            error_file.write(f"{header}\n\n")
            for error in errors:
                error_file.write(f"{error}\n")
            error_file.write("\n")


if __name__ == "__main__":
    main()
