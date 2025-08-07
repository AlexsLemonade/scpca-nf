#!/usr/bin/env python3

# Download reference files for the scpca-nf nextflow workflow to enable running
# the workflow without internet access by compute nodes. Optionally pulls
# container images for singularity or docker.
#
# Example usage:
# python3 get_refs.py --singularity


import argparse
import csv
import json
import os
import re
import shutil
import subprocess
import sys
import urllib.request
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--refdir",
        type=str,
        default="scpca-references",
        help="destination directory for downloaded reference files",
    )
    parser.add_argument(
        "--paramfile",
        type=str,
        default="localref_params.yaml",
        help="path to nextflow param file to write (default: `localref_params.yaml`)",
    )
    parser.add_argument(
        "--overwrite_refs",
        action="store_true",
        help="replace previously downloaded files",
    )
    parser.add_argument(
        "-r",
        "--revision",
        type=str,
        default="main",
        metavar="vX.X.X",
        help="tag for a specific workflow version (defaults to latest revision)",
    )
    parser.add_argument(
        "--star_index",
        action="store_true",
        help="get STAR index (required for genetic demultiplexing)",
    )
    parser.add_argument(
        "--cellranger_index",
        action="store_true",
        help="get Cell Ranger index (required for spatial data)",
    )
    parser.add_argument(
        "--docker", action="store_true", help="pull and cache images for docker"
    )
    parser.add_argument(
        "--singularity",
        action="store_true",
        help="pull and cache images for singularity",
    )
    parser.add_argument(
        "--singularity_dir",
        type=str,
        default="singularity",
        help="directory to store singularity image files (default: `singularity`)",
    )
    parser.add_argument(
        "--dry_run",
        action="store_true",
        help="perform a dry run; no large files will be downloaded",
    )
    args = parser.parse_args()

    # scpca-nf resource urls
    reffile_url = f"https://raw.githubusercontent.com/AlexsLemonade/scpca-nf/{args.revision}/config/reference_paths.config"
    refjson_url = f"https://raw.githubusercontent.com/AlexsLemonade/scpca-nf/{args.revision}/references/scpca-refs.json"

    containerfile_url = f"https://raw.githubusercontent.com/AlexsLemonade/scpca-nf/{args.revision}/config/containers.config"

    # download reference files
    print("Getting list of required reference files")
    try:
        ref_file = urllib.request.urlopen(reffile_url)
        json_file = urllib.request.urlopen(refjson_url)
    except urllib.error.URLError as e:
        print(e.reason)
        print("Reference file path download failed.")
        print(f"Is `{args.revision}` a valid release tag?")
        exit(1)

    # parse main reference file
    # gets all of the `param` variables that are set in `reference_paths.config`
    # and stores then in a dict
    ref_params = {}
    ref_re = re.compile(r'(?P<id>.+?)\s*=\s*([\'"])(?P<loc>.+)\2')
    for line in ref_file:
        match = ref_re.search(line.decode())
        if match:
            ref_params[match.group("id")] = match.group("loc")

    # get root location
    url_root = get_root_url(ref_params.get("ref_rootdir"))

    # create a list of all paths for all assemblies listed in reference json file
    json_paths = json.load(json_file)
    ref_paths = []
    for assembly_refs in json_paths.values():
        # add all the paths for single-file references
        ref_paths += get_single_files(assembly_refs)

        # get paths to salmon directories and add individual files
        ref_paths += get_salmon_files(assembly_refs)

        # if cellranger and star indices are requested, get individual files
        if args.cellranger_index:
            ref_paths += get_cellranger_files(assembly_refs)  # add cellranger files
        if args.star_index:
            ref_paths += get_star_files(assembly_refs)  # add star index files

    # barcode paths are still kept in the reference config file, so add those separately
    ref_paths += get_barcode_files(ref_params)

    # add paths to celltype reference files
    ref_paths += get_celltype_files(ref_params)

    # download all the files stored in `ref_paths`
    download_files(
        ref_paths,
        url_root,
        output_dir=args.refdir,
        overwrite=args.overwrite_refs,
        dry_run=args.dry_run,
    )

    # write param file if requested
    if args.paramfile and not args.dry_run:
        write_paramfile(args.paramfile, args.refdir, ref_params)

    # get containers if requested
    if args.docker or args.singularity:
        download_containers(
            containerfile_url,
            get_docker=args.docker,
            get_singularity=args.singularity,
            singularity_dir=args.singularity_dir,
            dry_run=args.dry_run,
        )
    print("All files downloaded successfully!")


def get_root_url(root_uri: str) -> str:
    # split out protocol from the root URI
    root_parts = root_uri.split("://", maxsplit=1)
    if root_parts[0] == "s3":
        # if S3, convert bucket path to https:// url
        bucket_path = root_parts[1].split("/", maxsplit=1)
        url_root = f"https://{bucket_path[0]}.s3.amazonaws.com"
        if len(bucket_path) > 1:
            url_root += f"/{bucket_path[1]}"
    elif root_parts[0] in ["http", "https", "ftp"]:
        # otherwise, just get the location
        url_root = root_uri
    else:
        raise ValueError(f"`{root_uri}` is not a supported remote location.")
    return url_root


def get_single_files(param_dict: dict[str, str]) -> list[Path]:
    # the keys here are the param variables that contain single files
    ref_keys = [
        "ref_fasta",
        "ref_fasta_index",
        "ref_gtf",
        "mito_file",
        "t2g_3col_path",
        "t2g_bulk_path",
    ]

    return [Path(param_dict.get(k)) for k in ref_keys]


def get_salmon_files(param_dict: dict[str, str]) -> list[Path]:
    # param variables that are salmon index directories
    salmon_keys = ["splici_index", "salmon_bulk_index"]

    # salmon index files within index dir (must be downloaded individually through http)
    salmon_index_files = [
        "complete_ref_lens.bin",
        "ctable.bin",
        "ctg_offsets.bin",
        "duplicate_clusters.tsv",
        "info.json",
        "mphf.bin",
        "pos.bin",
        "pre_indexing.log",
        "rank.bin",
        "ref_indexing.log",
        "refAccumLengths.bin",
        "reflengths.bin",
        "refseq.bin",
        "seq.bin",
        "versionInfo.json",
    ]

    salmon_paths = []
    sa_dir = [Path(param_dict.get(k)) for k in salmon_keys]
    for dir in sa_dir:
        salmon_paths += [dir / f for f in salmon_index_files]  # add salmon files

    return salmon_paths


def get_star_files(param_dict: dict[str, str]) -> list[Path]:
    # star index files within index dir (must be downloaded individually through http)
    star_index_files = [
        "chrLength.txt",
        "chrName.txt",
        "chrNameLength.txt",
        "chrStart.txt",
        "exonGeTrInfo.tab",
        "exonInfo.tab",
        "geneInfo.tab",
        "Genome",
        "genomeParameters.txt",
        "Log.out",
        "SA",
        "SAindex",
        "sjdbInfo.txt",
        "sjdbList.fromGTF.out.tab",
        "sjdbList.out.tab",
        "transcriptInfo.tab",
    ]
    star_dir = Path(param_dict.get("star_index"))
    return [star_dir / f for f in star_index_files]


def get_cellranger_files(param_dict: dict[str, str]) -> list[Path]:
    # Cell Ranger index files within index dir (must be downloaded individually through http)
    cr_index_files = [
        "reference.json",
        "fasta/genome.fa",
        "fasta/genome.fa.fai",
        "genes/genes.gtf.gz",
        "star/chrLength.txt",
        "star/chrName.txt",
        "star/chrNameLength.txt",
        "star/chrStart.txt",
        "star/exonGeTrInfo.tab",
        "star/exonInfo.tab",
        "star/geneInfo.tab",
        "star/Genome",
        "star/genomeParameters.txt",
        "star/SA",
        "star/SAindex",
        "star/sjdbInfo.txt",
        "star/sjdbList.fromGTF.out.tab",
        "star/sjdbList.out.tab",
        "star/transcriptInfo.tab",
    ]
    cr_dir = Path(param_dict.get("cellranger_index"))
    return [cr_dir / f for f in cr_index_files]


def get_barcode_files(param_dict: dict[str, str]) -> list[Path]:
    # barcode files on S3 within the barcode_dir (must be downloaded individually through http)
    barcode_files = [
        "3M-february-2018.txt",
        "737K-august-2016.txt",
        "cellranger_mit_license.txt",
        "visium-v1.txt",
        "visium-v2.txt",
    ]
    # regular expression for parameter expansion
    root_re = re.compile(r"^\$\{?(params.)?ref_rootdir\}?\/")
    # strip off root dir if it is there
    barcode_dir = root_re.sub("", param_dict.get("barcode_dir"))
    barcode_dir = Path(barcode_dir)
    return [barcode_dir / f for f in barcode_files]


def get_celltype_files(param_dict: dict[str, str]) -> list[Path]:
    # get celltype reference paths
    root_url = get_root_url(param_dict.get("ref_rootdir"))

    # regular expressions to match variables in parameter strings:
    # ${params.ref_rootdir}, $celltype_ref_dir, etc. ({} and `params.` prefix optional)
    root_re = re.compile(r"^\$\{?(params.)?ref_rootdir\}?")
    celltype_dir_re = re.compile(r"^\$\{?(params.)?celltype_ref_dir\}?")

    celltype_ref_dir = param_dict.get("celltype_ref_dir")
    # replace root placeholder if needed
    celltype_ref_url = root_re.sub(root_url, celltype_ref_dir)

    # get paths (with variables) to celltype reference files
    singler_model_dir = param_dict.get("singler_models_dir")
    cellassign_ref_dir = param_dict.get("cellassign_ref_dir")

    # replace celltype ref placeholders
    singler_model_dir = celltype_dir_re.sub(celltype_ref_url, singler_model_dir)
    cellassign_ref_dir = celltype_dir_re.sub(celltype_ref_url, cellassign_ref_dir)

    # alternatively, these might be tied to the root url
    singler_model_dir = root_re.sub(root_url, singler_model_dir)
    cellassign_ref_dir = root_re.sub(root_url, cellassign_ref_dir)

    # get celltype reference file tables from directories
    singler_model_url = f"{singler_model_dir}/singler_models.tsv"
    cellassign_ref_url = f"{cellassign_ref_dir}/cellassign_references.tsv"

    try:
        singler_model_file = urllib.request.urlopen(singler_model_url)
        cellassign_ref_file = urllib.request.urlopen(cellassign_ref_url)
    except urllib.error.URLError as e:
        print(e.reason)
        print("Celltype reference file path download failed.")
        raise

    singler_refs = csv.DictReader(
        (line.decode() for line in singler_model_file), delimiter="\t"
    )
    cellassign_refs = csv.DictReader(
        (line.decode() for line in cellassign_ref_file), delimiter="\t"
    )

    # get path strings to all celltype reference files
    celltype_ref_files = [
        f"{singler_model_dir}/{row['filename']}" for row in singler_refs
    ] + [f"{cellassign_ref_dir}/{row['filename']}" for row in cellassign_refs]

    # strip off root url and convert to path before returning
    celltype_ref_files = [
        Path(f.removeprefix(f"{root_url}/")) for f in celltype_ref_files
    ]
    return celltype_ref_files


def download_files(
    ref_paths: list, url_root: str, output_dir: Path, overwrite: bool, dry_run=False
) -> None:
    ## download all the files and put them in the correct locations ##
    print("Downloading reference files... (This might take a while)")
    for path in ref_paths:
        outfile = output_dir / path
        if outfile.exists() and not overwrite:
            continue
        print(f"Getting {path}")
        # make parents
        outfile.parent.mkdir(exist_ok=True, parents=True)
        # download and write
        file_url = f"{url_root}/{path}"
        if dry_run:
            print(f"Would download {file_url} to {outfile}")
        else:
            try:
                urllib.request.urlretrieve(file_url, outfile)
            except urllib.error.URLError as e:
                print(e.reason)
                print(
                    f"The file download failed for {file_url}; please check the URL for errors",
                    file=sys.stderr,
                )
                raise
    print(
        "Done with reference file downloads\n"
        f"Reference files can be found at '{output_dir}'\n"
    )


def write_paramfile(paramfile: str, refdir: str, param_dict: dict[str, str]) -> None:
    paramfile = Path(paramfile)
    # check if paramfile exists & move old if needed
    if paramfile.exists():
        print(
            f"A file already exists at `{paramfile}`; renaming previous file to `{paramfile.name}.bak`"
        )
        shutil.move(paramfile, str(paramfile) + ".bak")

    # always update the rootdir
    updated_params = {
        "ref_rootdir": os.path.abspath(refdir),
    }
    # Find any parameters that contain substitutions and fill in with substituted values
    nf_param_re = re.compile(r"\$\{?(params.)?(?P<variable>.+?)\}?/")
    for key, value in param_dict.items():
        match = nf_param_re.match(value)
        if match and match.group("variable") in updated_params:
            updated_params[key] = nf_param_re.sub(
                f"{updated_params[match.group('variable')]}/", value
            )
    # we don't want to change `singler_references_dir` as that is not downloaded
    updated_params.pop("singler_references_dir", None)

    with open(paramfile, "w") as f:
        f.write(
            "# Local Nextflow reference file parameters, generated by `get_refs.py`\n"
            "# To use with Nextflow, include the option `-params-file <path/to/this/file>` \n\n"
        )
        for key, value in updated_params.items():
            f.write(f"{key}: {value}\n")


def download_containers(
    containerfile_url: str,
    get_docker: bool,
    get_singularity: bool,
    singularity_dir: str,
    dry_run=False,
) -> None:
    ## Get docker containers from workflow
    print("Getting list of required containers")
    containers = {}
    try:
        container_file = urllib.request.urlopen(containerfile_url)
    except urllib.error.URLError as e:
        print(e.reason)
        print(
            f"The file download failed for {containerfile_url}; please check the URL for errors"
        )
        raise

    # pattern match to find container id & location
    container_re = re.compile(r'(?P<id>.+_CONTAINER)\s*=\s*([\'"])(?P<loc>.+)\2')
    for line in container_file:
        match = container_re.search(line.decode())
        if match:
            containers[match.group("id")] = match.group("loc")

    # pull docker images if requested
    if get_docker and not dry_run:
        print("Pulling docker images...")
        for loc in containers.values():
            subprocess.run(["docker", "pull", loc])
        print("Done pulling docker images\n")

    # pull singularity images if requested (to optionally specified cache location)
    if get_singularity and not dry_run:
        print("Pulling singularity images...")
        image_dir = Path(singularity_dir)
        image_dir.mkdir(parents=True, exist_ok=True)
        for loc in containers.values():
            # create image file name from location for nextflow
            image_file = loc.replace("/", "-").replace(":", "-") + ".img"
            print(image_file)
            image_path = image_dir / image_file
            print(image_path)
            if image_path.exists():
                image_path.unlink()
            subprocess.run(
                ["singularity", "pull", "--name", image_path, f"docker://{loc}"]
            )
        print("Done pulling singularity images")
        print(
            f"Singularity images located at {image_dir.absolute()}."
            " Be sure to set `singularity.cacheDir` in your configuration file to this value"
        )
        print()


if __name__ == "__main__":
    main()
