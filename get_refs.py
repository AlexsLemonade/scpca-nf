#!/usr/bin/env python3

import argparse
import os
import pathlib
import re
import shutil
import subprocess
import sys
import urllib.request

parser = argparse.ArgumentParser()
parser.add_argument("--refdir", type=str,
                    default="scpca-references",
                    help = "destination directory for downloaded reference files")
parser.add_argument("--replace",
                    action = "store_true",
                    help = "replace previously downloaded files")
parser.add_argument("--paramfile", type=str,
                    default="local_refs.yaml",
                    help = "nextflow param file to write (default: `local_refs.params`)")
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
parser.add_argument("--singularity_cache", type=str,
                    metavar = "CACHE_DIR",
                    help = "cache directory for singularity")
args = parser.parse_args()

# scpca-nf resource urls
aws_root = "https://scpca-references.s3.amazonaws.com"
containerfile_url = f"https://raw.githubusercontent.com/AlexsLemonade/scpca-nf/{args.revision}/config/containers.config"

# genome reference files
assembly = "Homo_sapiens.GRCh38.104"
genome_dir = pathlib.Path("homo_sapiens/ensembl-104")
ref_subdirs =[
    "fasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz",
    "fasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai",
    "annotation/Homo_sapiens.GRCh38.104.gtf.gz",
    "annotation/Homo_sapiens.GRCh38.104.mitogenes.txt",
    "annotation/Homo_sapiens.GRCh38.104.spliced_intron.tx2gene_3col.tsv",
    "annotation/Homo_sapiens.GRCh38.104.spliced_cdna.tx2gene.tsv"
]
ref_paths = [genome_dir / sd for sd in ref_subdirs]

# salmon index files
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

salmon_index_dirs = [
    genome_dir / "salmon_index/Homo_sapiens.GRCh38.104.spliced_intron.txome",
    genome_dir / "salmon_index/Homo_sapiens.GRCh38.104.spliced_cdna.txome"
]
for sa_dir in salmon_index_dirs:
    ref_paths += [sa_dir / f for f in salmon_index_files]

# star index files
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

star_dir = genome_dir / "star_index/Homo_sapiens.GRCh38.104.star_idx"
if args.star_index:
    ref_paths += [star_dir / f for f in star_index_files]

# Cell Ranger index files
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
cr_dir = genome_dir / "cellranger_index/Homo_sapiens.GRCh38.104_cellranger_full"
if args.cellranger_index:
    ref_paths += [cr_dir / f for f in cr_index_files]

# barcode file paths
barcode_dir = pathlib.Path("barcodes/10X")
barcode_files = [
    "3M-february-2018.txt",
    "737K-august-2016.txt",
    "cellranger_mit_license.txt",
    "visium-v1.txt",
    "visium-v2.txt"
]

ref_paths += [barcode_dir / f for f in barcode_files]

## download all the files and put them in the correct locations ##
print("Downloading reference files...")
for path in ref_paths[0:2]:
    outfile = args.refdir / path
    if outfile.exists() and not args.replace:
        continue
    print(f"Getting {path}")
    # make parents
    outfile.parent.mkdir(exist_ok=True, parents = True)
    # download and write
    file_url = f"{aws_root}/{path}"
    try:
        urllib.request.urlretrieve(file_url, outfile)
    except urllib.error.URLError as e:
        print(e.reason)
        print(f"The file download failed for {file_url}, please check the URL for errors",
              file = sys.stderr)
        exit(1)
print("Done with reference file downloads\n"
      f"Reference files can be found at '{args.refdir}'\n")

# write param file if requested
if args.paramfile:
    pfile = pathlib.Path(args.paramfile)
    # check if paramfile exists & move old if needed
    if pfile.exists():
        print(f"A file already exists at `{pfile}`, renaming previous file to `{pfile.name}.bak`")
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
        print(f"The file download failed for {container_url}, please check the URL for errors")
        print(f"Is `{args.revision}` a valid release tag?")
        exit(1)

    # pattern match to find container id & location
    container_re = re.compile(r'(?P<id>.+_CONTAINER)\s*=\s*([\'"])(?P<loc>.+)\2')
    for line in container_file:
        match = container_re.search(line.decode())
        if match:
            containers[match.group('id')] = match.group('loc')

# pull docker images ##
if args.docker:
    print("Pulling docker images...")
    for loc in containers.values():
        subprocess.run(["docker", "pull", loc])
    print("Done pulling docker images\n")

# pull singularity images (to optionally specified cache location)
if args.singularity:
    print("Pulling singularity images...")
    if args.singularity_cache:
        os.environ['SINGULARITY_CACHEDIR'] = os.path.abspath(args.singularity_cache)
    for loc in containers.values():
        subprocess.run(
            ["singularity", "pull", "--force", f"docker://{loc}"],
            env = os.environ
        )
    print("Done pulling singularity images")
    if args.singularity_cache:
        print(f"Singularity images located at {os.environ['SINGULARITY_CACHEDIR']}")
    print()
