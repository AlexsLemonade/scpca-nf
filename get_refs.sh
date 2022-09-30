#!/bin/bash
set -euo pipefail

dest_dir=${1:-"scpca-references"}


containerfile_url="https://raw.githubusercontent.com/AlexsLemonade/scpca-nf/main/config/containers.config"
aws_root="https://scpca-references.s3.amazonaws.com"

ref_dir="homo_sapiens/ensembl-104"
ref_paths=(
    "${ref_dir}/fasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
    "${ref_dir}/fasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai"
    "${ref_dir}/annotation/Homo_sapiens.GRCh38.104.gtf.gz"
    "${ref_dir}/annotation/Homo_sapiens.GRCh38.104.mitogenes.txt"
    "${ref_dir}/annotation/Homo_sapiens.GRCh38.104.spliced_intron.tx2gene_3col.tsv"
    "${ref_dir}/annotation/Homo_sapiens.GRCh38.104.spliced_cdna.tx2gene.tsv"
)

salmon_index_files=(
    "complete_ref_lens.bin"
    "ctable.bin"
    "ctg_offsets.bin"
    "duplicate_clusters.tsv"
    "info.json"
    "mphf.bin"
    "pos.bin"
    "pre_indexing.log"
    "rank.bin"
    "ref_indexing.log"
    "refAccumLengths.bin"
    "reflengths.bin"
    "refseq.bin"
    "seq.bin"
    "versionInfo.json"
)

salmon_index_dirs=(
    "${ref_dir}/salmon_index/Homo_sapiens.GRCh38.104.spliced_intron.txome"
    "${ref_dir}/salmon_index/Homo_sapiens.GRCh38.104.spliced_cdna.txome"
)
for dir in ${salmon_index_dirs[@]}
do
    for file in ${salmon_index_files[@]}
    do
        ref_paths+=("${dir}/${file}")
    done
done


star_index_files=(
    "chrLength.txt"
    "chrName.txt"
    "chrNameLength.txt"
    "chrStart.txt"
    "exonGeTrInfo.tab"
    "exonInfo.tab"
    "geneInfo.tab"
    "Genome"
    "genomeParameters.txt"
    "Log.out"
    "SA"
    "SAindex"
    "sjdbInfo.txt"
    "sjdbList.fromGTF.out.tab"
    "sjdbList.out.tab"
    "transcriptInfo.tab"
)

star_dir="${ref_dir}/star_index/Homo_sapiens.GRCh38.104.star_idx"
for file in ${star_index_files[@]}
do
    ref_paths+=("${star_dir}/${file}")
done


barcode_files=(
    "3M-february-2018.txt"
    "737K-august-2016.txt"
    "cellranger_mit_license.txt"
    "visium-v1.txt"
    "visium-v1.txt"
)

barcode_dir="barcodes/10X"
for file in ${barcode_files[@]}
do
    ref_paths+=("${barcode_dir}/${file}")
done

for path in ${ref_paths[@]}
do
    echo "Getting $path"
    curl -s --create-dirs "$aws_root/$path" -o "$dest_dir/$path"
done


containers=`curl -s https://raw.githubusercontent.com/AlexsLemonade/scpca-nf/main/config/containers.config \
    | grep CONTAINER \
    | cut -d"'" -f 2 \
    | grep -v "^$" `
