#!/usr/bin/env python3

import argparse
import pathlib
from signal import SIG_DFL
import sys
import urllib.request

parser = argparse.ArgumentParser()
parser.add_argument("--refdir", type=str,
                    default="scpca-references",
                    help = "destination directory for downloaded reference files")
parser.add_argument("--replace", action = argparse.BooleanOptionalAction,
                    default = False,
                    help = "replace previously downloaded files")
parser.add_argument("--revision", type=str,
                    default="main",
                    metavar = "vX.X.X",
                    help = "tag for a specific workflow version (defaults to latest revision)")

args = parser.parse_args()

aws_root = "https://scpca-references.s3.amazonaws.com"
containerfile_url = f"https://raw.githubusercontent.com/AlexsLemonade/scpca-nf/{args.revision}/config/containers.config"

# genome reference files
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

star_dir = genome_dir/ "star_index/Homo_sapiens.GRCh38.104.star_idx"
ref_paths += [star_dir / f for f in star_index_files]

# get barcode file paths
barcode_dir = pathlib.Path("barcodes/10X")
barcode_files = [
    "3M-february-2018.txt",
    "737K-august-2016.txt",
    "cellranger_mit_license.txt",
    "visium-v1.txt",
    "visium-v2.txt"
]

ref_paths += [barcode_dir / f for f in barcode_files]

# download all the files and put them in the correct locations
print("Downloading reference files:")
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
    except urllib.error.URLError:
        print(f"The file download failed for {file_url}, please check the URL for errors",
              file = sys.stderr)
        exit(1)

