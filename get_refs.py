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
    (
        parser.add_argument(
            "--singularity_dir",
            type=str,
            default="singularity",
            help="directory to store singularity image files (default: `singularity`)",
        ),
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
    celltype_metadata_url = f"https://raw.githubusercontent.com/AlexsLemonade/scpca-nf/{args.revision}/references/celltype-reference-metadata.tsv"

    containerfile_url = f"https://raw.githubusercontent.com/AlexsLemonade/scpca-nf/{args.revision}/config/containers.config"

    # download reference files
    print("Getting list of required reference files")
    try:
        ref_file = urllib.request.urlopen(reffile_url)
        json_file = urllib.request.urlopen(refjson_url)
        celltype_metadata_file = urllib.request.urlopen(celltype_metadata_url)
    except urllib.error.URLError as e:
        print(e.reason)
        print("Reference file path download failed.")
        print(f"Is `{args.revision}` a valid release tag?")
        exit(1)

    # parse main reference file
    # gets all of the `param` variables that are set in `reference_paths.config`
    # and stores then in a dict
    refs = {}
    ref_re = re.compile(r'(?P<id>.+?)\s*=\s*([\'"])(?P<loc>.+)\2')
    for line in ref_file:
        match = ref_re.search(line.decode())
        if match:
            refs[match.group("id")] = match.group("loc")

    # get root location
    url_root = get_root_url(refs.get("ref_rootdir"))

    # create a list of all paths for all assemblies listed in reference json file
    json_paths = json.load(json_file)
    ref_paths = []
    for assembly_refs in json_paths.values():
        # add all the paths for single-file references
        ref_paths += get_single_files(assembly_refs)

        # get paths to salmon directories and add individual files
        ref_paths += get_salmon_files(assembly_refs)

        # if cellranger and star indicies are requested, get individual files
        if args.cellranger_index:
            ref_paths += get_cellranger_files(assembly_refs)  # add cellranger files
        if args.star_index:
            ref_paths += get_star_files(assembly_refs)  # add star index files

    # barcode paths are still kept in the reference config file, so add those separately
    ref_paths += get_barcode_files(refs)

    # paths to celltype reference files
    # read in celltype reference metadata
    celltype_metadata = csv.DictReader(
        (line.decode() for line in celltype_metadata_file), delimiter="\t"
    )
    ref_paths += get_celltype_files(refs, celltype_metadata)

    download_files(
        ref_paths,
        url_root,
        output_dir=args.refdir,
        overwrite=args.overwrite_refs,
        dry_run=args.dry_run,
    )

    # write param file if requested
    if args.paramfile and not args.dry_run:
        pfile = Path(args.paramfile)
        # check if paramfile exists & move old if needed
        if pfile.exists():
            print(
                f"A file already exists at `{pfile}`; renaming previous file to `{pfile.name}.bak`"
            )
            shutil.move(pfile, str(pfile) + ".bak")
        # create parameter dictionary
        nf_params = {
            "ref_rootdir": os.path.abspath(args.refdir),
        }
        with open(pfile, "w") as f:
            f.write(
                "# local nextflow reference file parameters, generated by `get_refs.py`\n\n"
            )
            for key, value in nf_params.items():
                f.write(f"{key}: {value}\n")

    ## Get docker containers from workflow
    if args.singularity or args.docker:
        print("Getting list of required containers")
        containers = {}
        try:
            container_file = urllib.request.urlopen(containerfile_url)
        except urllib.error.URLError as e:
            print(e.reason)
            print(
                f"The file download failed for {containerfile_url}; please check the URL for errors"
            )
            print(f"Is `{args.revision}` a valid release tag?")
            exit(1)

        # pattern match to find container id & location
        container_re = re.compile(r'(?P<id>.+_CONTAINER)\s*=\s*([\'"])(?P<loc>.+)\2')
        for line in container_file:
            match = container_re.search(line.decode())
            if match:
                containers[match.group("id")] = match.group("loc")

    # pull docker images if requested
    if args.docker and not args.dry_run:
        print("Pulling docker images...")
        for loc in containers.values():
            subprocess.run(["docker", "pull", loc])
        print("Done pulling docker images\n")

    # pull singularity images if requested (to optionally specified cache location)
    if args.singularity and not args.dry_run:
        print("Pulling singularity images...")
        image_dir = Path(args.singularity_dir)
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
        raise ValueError("`ref_rootdir` is not a supported remote location.")
    return url_root


def get_single_files(ref_dict: dict[str, str]) -> list[Path]:
    # the keys here are the param variables that contain single files
    ref_keys = [
        "ref_fasta",
        "ref_fasta_index",
        "ref_gtf",
        "mito_file",
        "t2g_3col_path",
        "t2g_bulk_path",
    ]

    return [Path(ref_dict.get(k)) for k in ref_keys]


def get_salmon_files(ref_dict: dict[str, str]) -> list[Path]:
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
    sa_dir = [Path(ref_dict.get(k)) for k in salmon_keys]
    for dir in sa_dir:
        salmon_paths += [dir / f for f in salmon_index_files]  # add salmon files

    return salmon_paths


def get_star_files(ref_dict: dict[str, str]) -> list[Path]:
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
    star_dir = Path(ref_dict.get("star_index"))
    return [star_dir / f for f in star_index_files]


def get_cellranger_files(ref_dict: dict[str, str]) -> list[Path]:
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
    cr_dir = Path(ref_dict.get("cellranger_index"))
    return [cr_dir / f for f in cr_index_files]


def get_barcode_files(ref_dict: dict[str, str]) -> list[Path]:
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
    barcode_dir = root_re.sub("", ref_dict.get("barcode_dir"))
    barcode_dir = Path(barcode_dir)
    return [barcode_dir / f for f in barcode_files]


def get_celltype_files(
    ref_dict: dict[str, str], celltype_metadata: dict[str, str]
) -> list[Path]:
    # get celltype reference paths
    root_url = get_root_url(ref_dict.get("ref_rootdir"))

    root_re = re.compile(r"^\$\{?(params.)?ref_rootdir\}?")
    celltype_dir_re = re.compile(r"^\$\{?(params.)?celltype_ref_dir\}?")

    celltype_ref_dir = ref_dict.get("celltype_ref_dir")
    # replace root placeholder if needed
    celltype_ref_url = root_re.sub(root_url, celltype_ref_dir)

    # get paths (with variables) to celltype reference files
    singler_ref_dir = ref_dict.get("singler_references_dir")
    singler_model_dir = ref_dict.get("singler_models_dir")
    cellassign_ref_dir = ref_dict.get("cellassign_ref_dir")

    # replace celltype ref placeholders
    singler_ref_dir = celltype_dir_re.sub(celltype_ref_url, singler_ref_dir)
    singler_model_dir = celltype_dir_re.sub(celltype_ref_url, singler_model_dir)
    cellassign_ref_dir = celltype_dir_re.sub(celltype_ref_url, cellassign_ref_dir)

    # alternatively, these might be tied to the root url
    singler_ref_dir = root_re.sub(root_url, singler_ref_dir)
    singler_model_dir = root_re.sub(root_url, singler_model_dir)
    cellassign_ref_dir = root_re.sub(root_url, cellassign_ref_dir)

    # get versions of celltype reference files
    celldex_version_url = f"{singler_ref_dir}/celldex_version.txt"
    panglao_version_url = f"{cellassign_ref_dir}/PanglaoDB_version.txt"

    celldex_version = (
        urllib.request.urlopen(celldex_version_url)
        .read()
        .decode()
        .strip()
        .replace(".", "-")  # replace periods with dashes for file names
    )
    panglao_version = (
        urllib.request.urlopen(panglao_version_url).read().decode().strip()
    )

    celltype_ref_files = []
    for row in celltype_metadata:
        ref_name = row["celltype_ref_name"]
        method = row["celltype_method"]
        source = row["celltype_ref_source"]
        if method.lower() == "singler" and source == "celldex":
            celltype_ref_files.append(
                f"{singler_model_dir}/{ref_name}_{source}_{celldex_version}_model.rds"
            )
        elif method.lower() == "cellassign" and source == "PanglaoDB":
            celltype_ref_files.append(
                f"{cellassign_ref_dir}/{ref_name}_{source}_{panglao_version}.tsv"
            )
        else:
            raise ValueError(
                f"Celltype reference file for {method} and {source} is not supported."
            )
    # strip off root url and convert to path before returning
    celltype_ref_files = [
        Path(f.removeprefix(f"{root_url}/")) for f in celltype_ref_files
    ]
    return celltype_ref_files


def download_files(
    ref_paths: list, url_root: str, output_dir: Path, overwrite: bool, dry_run=False
):
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


if __name__ == "__main__":
    main()
