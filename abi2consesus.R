#! /usr/bin/env Rscript

require(ape)
require(reshape2)
require(phangorn)
require(stringi)
require(stringr)
require(argparse)
require(ggplot2)
require(here)
require(grid)
require(gridExtra)
# Bioconductor packages
require(DECIPHER)
require(Biostrings)
require(sangerseqR)
#github("roblanf/sangeranalyseR")
require(sangeranalyseR)

args <- parser$parse_args()
input_folder <-args$input_folder
forward_suffix <- args$forward_suffix
reverse_suffix <- args$reverse_suffix

parser <- ArgumentParser()
parser$add_argument("input_folder", help="input folder where abi files are located.")
parser$add_argument("-f","--foward-suffix", default="_F.ab1", dest= "forward_suffix")
parser$add_argument("-r","--reverse-suffix", default="_R.ab1", dest= "reverse_suffix")

consesus = make.consensus.seqs(input_folder,
                               forward_suffix, 
                               reverse_suffix,
                               trim = TRUE, 
                               trim.cutoff = 0.015)


consesus$read.summaries
consesus$consensus.summaries
consesus$consensus.sequences
BrowseSeqs(consesus$consensus.alignment)
plot(cs$consensus.tree)
fname = paste0(basename(input_folder), "_consesus",".fasta")
write.dna(cs$consensus.sequences, 
          file = fname,
          format = 'fasta',
          nbcol = -1, 
          colsep = "",
          colw = 10000000)
