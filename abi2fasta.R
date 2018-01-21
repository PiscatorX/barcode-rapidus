#! /usr/bin/env Rscript

require(argparse)
cwd = getwd()

parser <- ArgumentParser(description="Process sanger ab1 reads and generate fasta files")
parser$add_argument("input_folder", help="input folder where abi files are located")
parser$add_argument("-f","--foward-suffix", default="_F.ab1", dest= "forward_suffix")
parser$add_argument("-r","--reverse-suffix", default="_R.ab1", dest= "reverse_suffix")
parser$add_argument("-o","--output-folder", dest="output_folder", default=cwd,
                    help="directory to write output files")
parser$add_argument("-s","--split-char", dest="split_char", default="_",
                    help="character to split on to obtain sequence id from filename eg. for filname: 'LPA1.1_RH1-rbcL_1_R.ab1' spliton '_' (default)")
parser$add_argument("-t","--trim-cutoff", dest="trim_cutoff", default=0.015,
                    help="cutoff for sequencing trimming")
args <- parser$parse_args()

require(ape)
require(reshape2)
require(phangorn)
require(stringi)
require(stringr)
require(ggplot2)
#require(here)
require(grid)
require(gridExtra)
# Bioconductor packages
require(DECIPHER)
require(Biostrings)
require(sangerseqR)
require(tools)
#github("roblanf/sangeranalyseR")
require(sangeranalyseR)




setwd(cwd)
input_folder <- normalizePath(args$input_folder)
input_folder <- "/home/drewx/Documents/sanpccX/TestDev"
forward_suffix <- args$forward_suffix
reverse_suffix <- args$reverse_suffix
split_char <- args$split_char


read_sets <- make.readsets(input_folder, 
              forward_suffix, 
              reverse_suffix, 
              trim = TRUE, 
              trim.cutoff = 0.015)


base_seq_ids = character(length(read_sets$readsets))
i = 1
for (seq_name in names(read_sets$readsets)){
  seq_id <- basename(seq_name)
  new_seq_id <- paste(strsplit(seq_id, split_char)[[1]][[1]], seq_id)
  base_seq_ids[i] <- new_seq_id
  i <- i + 1
}

names(read_sets$readsets) <- base_seq_ids

read_sets$readsets
read_sets$read.summaries

fname = paste0(basename(input_folder),"_reads",".fasta")
setwd(args$output_folder)
write.dna(read_sets$readsets,
          file = fname,
          format = 'fasta', 
          nbcol = -1, 
          colsep = "", 
          colw = 10000000)

