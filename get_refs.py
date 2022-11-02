#!/usr/bin/env python3

# Download reference files for the scpca-nf nextflow workflow to enable running
# the workflow without internet access by compute nodes. Optionally pulls
# container images for singularity or docker.
#
# Example usage:
# python3 get_refs.py --singularity


import argparse
import os
import re
import shutil
import subprocess
import sys
import urllib.request

from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--refdir", type=str,
                    default="scpca-references",
                    help = "destination directory for downloaded reference files")
parser.add_argument("--paramfile", type=str,
                    default="localref_params.yaml",
                    help = "path to nextflow param file to write (default: `localref_params.yaml`)")
parser.add_argument("--overwrite_refs",
                    action = "store_true",
                    help = "replace previously downloaded files")
parser.add_argument("--revision", type=str,
                    default="main",
                    metavar = "vX.X.X",
                    help = "tag for a specific workflow version (defaults to latest revision)")
parser.add_argument("--star_index",
                    action = "store_true",
                    help = "get STAR index (required for genetic demultiplexing)")
parser.add_argument("--cellranger_index",
                    action = "store_true",
                    help = "get Cell Ranger index (required for spatial data)")
parser.add_argument("--docker",
                    action = "store_true",
                    help = "pull and cache images for docker")
parser.add_argument("--singularity",
                    action = "store_true",
                    help = "pull and cache images for singularity")
parser.add_argument("--singularity_dir", type=str,
                    default = "singularity",
                    help = "directory to store singularity image files (default: `singularity`)")
args = parser.parse_args()

# scpca-nf resource urls
reffile_url = f"https://raw.githubusercontent.com/AlexsLemonade/scpca-nf/{args.revision}/config/reference_paths.config"
containerfile_url = f"https://raw.githubusercontent.com/AlexsLemonade/scpca-nf/{args.revision}/config/containers.config"

# download reference file
print("Getting list of required reference files")
try:
    ref_file =  urllib.request.urlopen(reffile_url)
except urllib.error.URLError as e:
    print(e.reason)
    print(f"The file download failed for {reffile_url}; please check the URL for errors")
    print(f"Is `{args.revision}` a valid release tag?")
    exit(1)

# parse reference file
# gets all of the `param` variables that are set in `reference_paths.config`
# and stores then in a dict
refs = {}
ref_re = re.compile(r'(?P<id>.+?)\s*=\s*([\'"])(?P<loc>.+)\2')
for line in ref_file:
    match = ref_re.search(line.decode())
    if match:
        refs[match.group('id')] = match.group('loc')

# regular expressions for parameter expansion
root_re = re.compile(r'\$\{?(params.)?ref_rootdir\}?$')
refdir_re = re.compile(r'\$\{?(params.)?ref_dir\}?$')

# get assembly and root location
assembly = refs.get("assembly", "NA")
# split out protocol from the root URI
root_parts = refs.get("ref_rootdir").split('://', maxsplit = 1)
if root_parts[0] == 's3':
    # if S3, convert bucket path to https:// url
    bucket_path = root_parts[1].split("/", maxsplit = 1)
    url_root = f"https://{bucket_path[0]}.s3.amazonaws.com"
    if len(bucket_path) > 1:
        url_root += f"/{bucket_path[1]}"
elif root_parts[0] in ['http', 'https', 'ftp']:
    # otherwise, just get the location
    url_root = refs.get("ref_rootdir")
else:
    print("`ref_rootdir` is not a supported remote location.")
    exit(1)



# set the base directory (usually corresponding to a genome version)
genome_dir = Path(refs.get("ref_dir"))
# remove the first element if it is a variable
if root_re.match(genome_dir.parts[0]):
    genome_dir = genome_dir.relative_to(genome_dir.parts[0])

# single-file references
# the keys here are the param variables we will be downloading
ref_keys =[
    "ref_fasta",
    "ref_fasta_index",
    "ref_gtf",
    "mito_file",
    "t2g_3col_path",
    "t2g_bulk_path"
]
ref_paths = [Path(refs.get(k)) for k in ref_keys]
# replace initial part of path if it is `$params.ref_dir` or similar
ref_paths = [genome_dir / p.relative_to(p.parts[0])
             if refdir_re.match(p.parts[0]) else p
             for p in ref_paths]

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
    "versionInfo.json"
]
# param variables that are salmon index directories
salmon_keys = [
    "splici_index",
    "bulk_index"
]
for k in salmon_keys:
    sa_dir = Path(refs.get(k))
    if refdir_re.match(sa_dir.parts[0]):
        sa_dir = genome_dir / sa_dir.relative_to(sa_dir.parts[0])
    ref_paths += [sa_dir / f for f in salmon_index_files]

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
    "transcriptInfo.tab"
]
star_dir = Path(refs.get("star_index"))
if refdir_re.match(star_dir.parts[0]):
    star_dir = genome_dir / star_dir.relative_to(star_dir.parts[0])
if args.star_index:
    ref_paths += [star_dir / f for f in star_index_files]

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
    "star/transcriptInfo.tab"
]
cr_dir = Path(refs.get("cellranger_index"))
if refdir_re.match(cr_dir.parts[0]):
    cr_dir = genome_dir / cr_dir.relative_to(cr_dir.parts[0])
if args.cellranger_index:
    ref_paths += [cr_dir / f for f in cr_index_files]

# barcode files on S3 within the barcode_dir (must be downloaded individually through http)
barcode_files = [
    "3M-february-2018.txt",
    "737K-august-2016.txt",
    "cellranger_mit_license.txt",
    "visium-v1.txt",
    "visium-v2.txt"
]
barcode_dir = Path(refs.get("barcode_dir"))
if root_re.match(barcode_dir.parts[0]):
    barcode_dir = barcode_dir.relative_to(barcode_dir.parts[0])
ref_paths += [barcode_dir / f for f in barcode_files]

## download all the files and put them in the correct locations ##
print("Downloading reference files... (This might take a while)")
for path in ref_paths:
    outfile = args.refdir / path
    if outfile.exists() and not args.overwrite_refs:
        continue
    print(f"Getting {path}")
    # make parents
    outfile.parent.mkdir(exist_ok=True, parents = True)
    # download and write
    file_url = f"{url_root}/{path}"
    try:
        urllib.request.urlretrieve(file_url, outfile)
    except urllib.error.URLError as e:
        print(e.reason)
        print(f"The file download failed for {file_url}; please check the URL for errors",
              file = sys.stderr)
        exit(1)
print("Done with reference file downloads\n"
      f"Reference files can be found at '{args.refdir}'\n")

# write param file if requested
if args.paramfile:
    pfile = Path(args.paramfile)
    # check if paramfile exists & move old if needed
    if pfile.exists():
        print(f"A file already exists at `{pfile}`; renaming previous file to `{pfile.name}.bak`")
        shutil.move(pfile, str(pfile) + ".bak")
    # create parameter dictionary
    nf_params = {
        'assembly': assembly,
        'ref_rootdir': os.path.abspath(args.refdir)
    }
    with open(pfile, 'w') as f:
        f.write("# local nextflow reference file parameters, generated by `get_refs.py`\n\n")
        for key, value in nf_params.items():
            f.write(f"{key}: {value}\n")

## Get docker containers from workflow
if args.singularity or args.docker:
    print("Getting list of required containers")
    containers = {}
    try:
        container_file =  urllib.request.urlopen(containerfile_url)
    except urllib.error.URLError as e:
        print(e.reason)
        print(f"The file download failed for {container_url}; please check the URL for errors")
        print(f"Is `{args.revision}` a valid release tag?")
        exit(1)

    # pattern match to find container id & location
    container_re = re.compile(r'(?P<id>.+_CONTAINER)\s*=\s*([\'"])(?P<loc>.+)\2')
    for line in container_file:
        match = container_re.search(line.decode())
        if match:
            containers[match.group('id')] = match.group('loc')

# pull docker images if requested
if args.docker:
    print("Pulling docker images...")
    for loc in containers.values():
        subprocess.run(["docker", "pull", loc])
    print("Done pulling docker images\n")

# pull singularity images if requested (to optionally specified cache location)
if args.singularity:
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
        subprocess.run([
            "singularity", "pull",
            "--name", image_path,
            f"docker://{loc}"
        ])
    print("Done pulling singularity images")
    print(f"Singularity images located at {image_dir.absolute()}")
    print()
