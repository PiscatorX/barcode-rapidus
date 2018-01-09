#! /usr/bin/env Rscript

require(argparse)
cwd = getwd()
parser <- ArgumentParser()
parser <- ArgumentParser(description="Process sanger ab1 reads and generate fasta files")
parser$add_argument("input_folder", help="input folder where abi files are located.")
parser$add_argument("-g","--glob", dest="wild_card", help="sequence read glob wildcard")
parser$add_argument("-o","--output-folder", dest="output_folder", default=cwd,
                    help="directory to write output files")

args <- parser$parse_args()

require(ape)
require(reshape2)
require(phangorn)
require(stringi)
require(stringr)
require(ggplot2)
require(here)
require('grid')
require(gridExtra)
# Bioconductor packages
require(DECIPHER)
require(Biostrings)
require(sangerseqR)
# github("roblanf/sangeranalyseR")
require(sangeranalyseR)


input_folder =  args$input_folder
wild_card = args$wild_card

path_glob = file.path(input_folder, wild_card)
ab1_files = Sys.glob(path_glob)
output_folder = args$output_folder
reads_set = DNAStringSet()

print(output_folder)
quit()

output_folder = "/home/drewx/Documents/sanpccX/DevOps/AlignDevX"
for (ab1 in ab1_files){
    seq_abif = read.abif(ab1)
    seq_sanger = sangerseq(seq_abif)
    trims = trim.mott(seq_abif, cutoff=0.015)
    seq_untrimmed = seq_abif@data$PBAS.2 
    seq_trimmed = substring(seq_untrimmed, trims$start, trims$finish)
    print(paste("Untrimmed Sequence length:", nchar(seq_abif@data$PBAS.2)))
    print(paste("Trimmed Sequence length:", nchar(seq_trimmed))) 
    reads_set = append(reads_set, DNAStringSet(seq_trimmed)) 
    sp = secondary.peaks(seq_sanger)
    sp$secondary.peaks
    
}

fname_prefix = strsplit(basename(ab1_files)[1], '_')[[1]][1]
consesus_fname = paste0(file.prefix=fname_prefix, '_consesus','.fasta')
names(reads_set) = c(ab1_files)
merged_reads_set = merge.reads(reads_set)
names(merged_reads_set)
merged_reads_set
consesus_seq <- merged_reads_set$consensus
BrowseSeqs(merged_reads_set$alignment)
merged_reads_set$secondary.peak.columns
merged_reads_set$differeces
setwd(args$output_folder)
writeXStringSet(DNAStringSet(consesus_seq), consesus_fname, format="fasta")



