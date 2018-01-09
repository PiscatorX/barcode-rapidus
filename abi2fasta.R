#! /usr/bin/env Rscript


cwd = getwd()
require(argparse)
parser <- ArgumentParser(description="Process sanger ab1 reads and generate fasta files")
parser$add_argument("input_folder", help="input folder where abi files are located.")
parser$add_argument("-f","--foward-suffix", default="_F.ab1", dest= "forward_suffix")
parser$add_argument("-r","--reverse-suffix", default="_R.ab1", dest= "reverse_suffix")
parser$add_argument("-o","--output-folder", dest="output_folder", default=cwd,
                    help="directory to write output files")
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
input_folder <- args$input_folder
forward_suffix <- args$forward_suffix
reverse_suffix <- args$reverse_suffix

read_sets <- make.readsets(input_folder, 
              forward_suffix, 
              reverse_suffix, 
              trim = TRUE, 
              trim.cutoff = 0.015)

read_sets$readsets
read_sets$read.summaries

fname = paste0(basename(input_folder), "_reads",".fasta")
setwd(args$output_folder)
write.dna(read_sets$readsets,
          file = fname,
          format = 'fasta', 
          nbcol = -1, 
          colsep = "", 
          colw = 10000000)


