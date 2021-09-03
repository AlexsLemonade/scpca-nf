#!/usr/bin/env Rscript

# This script takes the input reference fasta of the genome and the 
# corresponding gtf file, identifies the regions of interest corresponding
# to spliced transcripts and introns, then subsets the genome and gtf for 
# only those particular regions, and finally outputs a spliced_intron.txome.fa, 
# spliced_intron.txome.gtf, and spliced_intron.txome.tx2gene.tsv all corresponding to the genomic 
# regions with spliced transcripts and introns

# This script is also used to generate the 3 column spliced_intron.tx2gene_3col.tsv needed for alevin-fry.

# finally, a list of mitochondrial genes is output for QC

# load needed packages
library(Biostrings)
library(BSgenome)
library(eisaR)
library(GenomicFeatures)
library(magrittr)
library(optparse)
library(tidyverse)

# Set up optparse options
option_list <- list(
  make_option(
    opt_str = c("-f", "--gtf"),
    type = "character",
    default = "annotation/Homo_sapiens.GRCh38.104.gtf.gz",
    help = "File path for input gtf file",
  ),
  make_option(
    opt_str = c("-g", "--genome"),
    type = "character",
    default = "fasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz",
    help = "File path for reference fasta file",
  ), 
  make_option(
    opt_str = c("-a", "--annotation_output"),
    type = "character",
    default = "annotation",
    help = "Directory to write output gtf, tx2gene, and mitochondrial gene list files",
  ), 
  make_option(
    opt_str = c("-o", "--fasta_output"),
    type = "character",
    default = "fasta",
    help = "Directory to write output transcript fasta files",
  ),
  make_option(
    opt_str = c("-s", "--assembly"),
    type = "character",
    default = ""
    help = "Prefix name containing organism and ensembl assembly version to be used for file naming"
  )
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

## files with spliced + intron regions
spliced_intron_fasta_file <- paste0(opt$assembly, ".spliced_intron.txome.fa.gz")
spliced_intron_gtf_file <- paste0(opt$assembly, ".spliced_intron.txome.gtf")
spliced_intron_tx2gene_file <- paste0(opt$assembly, ".spliced_intron.tx2gene.tsv")
spliced_intron_tx2gene_3col_file <- paste0(opt$assembly, ".spliced_intron.tx2gene_3col.tsv")

mito_file <- paste0(opt$assembly, ".mitogenes.txt")

# make final output file names needed

spliced_intron_fasta <- file.path(opt$fasta_output, spliced_intron_fasta_file)
spliced_intron_gtf <- file.path(opt$annotation_output, spliced_intron_gtf_file)
spliced_intron_tx2gene <- file.path(opt$annotation_output, 
                                    spliced_intron_tx2gene_file)

spliced_intron_tx2gene_3col <- file.path(opt$annotation_output, 
                                 spliced_intron_tx2gene_3col_file)

mito_out <- file.path(opt$annotation_output, mito_file)

# Check for output directory 
if (!dir.exists(opt$fasta_output)) {
  dir.create(opt$fasta_output)
}

if (!dir.exists(opt$annotation_output)) {
  dir.create(opt$annotation_output)
}


# extract GRanges object containing genomic coordinates of each annotated transcript
# both spliced transcripts and intronic regions
grl <- eisaR::getFeatureRanges(
  gtf = opt$gtf,
  featureType = c("spliced", "intron"), 
  flankLength = 86L,
  joinOverlappingIntrons = FALSE,
  verbose = TRUE
)

# load in primary genome
genome <- Biostrings::readDNAStringSet(opt$genome)
names(genome) <- stringr::word(names(genome), 1)


# add the levels and lengths
seqlevels(grl) <- seqlevels(genome)
seqlengths(grl) <- seqlengths(genome)

# trim grl used for splici
splici_grl <- trim(grl)

# extract seqs from trimmed grl for splici index
splici_seqs <- GenomicFeatures::extractTranscriptSeqs(
  x = genome, 
  transcripts = splici_grl
)

# remove sequence duplicates
splici_seqs <- unique(splici_seqs)
splici_grl <- splici_grl[names(splici_seqs)]

# write splici to fasta 
Biostrings::writeXStringSet(
  splici_seqs, filepath = spliced_intron_fasta, compress = TRUE
)

# write the associated annotations to gtf 
eisaR::exportToGtf(
  splici_grl, 
  filepath = spliced_intron_gtf
)

# create text file mapping transcript and intron identifiers to corresponding gene identifiers
# make 2 column Tx2Gene for all spliced and intron sequences
splici_tx2gene_df <- eisaR::getTx2Gene(splici_grl)

# write out 2 column Tx2 gene mapping
readr::write_tsv(splici_tx2gene_df, spliced_intron_tx2gene, col_names = FALSE)

# make 3 column Tx2 gene needed for alevin-fry USA mode
splici_tx2gene_df_3col <- splici_tx2gene_df  %>%
  # add status column
  # remove extra -I* on gene ID
  mutate(gene_id = stringr::word(gene_id,1,sep = '-'),
         #if transcript_id contains an added -I*, then it is usnpliced
         status = ifelse(stringr::str_detect(transcript_id, '-I'), 'U', 'S'))

# write 3 column tx2gene
readr::write_tsv(splici_tx2gene_df_3col, spliced_intron_tx2gene_3col, col_names = FALSE)


# reimport gtf to get list of mito genes 
gtf <- rtracklayer::import(spliced_intron_gtf)
mitogenes <- gtf[seqnames(gtf) == 'MT']

# write out mitochondrial gene list
writeLines(mitogenes$gene_id, mito_out)
