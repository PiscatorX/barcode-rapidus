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

#github("roblanf/sangeranalyseR")
require(sangeranalyseR)


setwd(cwd)
input_folder <- normalizePath(args$input_folder)
forward_suffix <- args$forward_suffix
reverse_suffix <- args$reverse_suffix
split_char <- args$split_char

read_sets <- make.consensus.seqs(input_folder, 
              forward_suffix, 
              reverse_suffix,
              min.reads = 0,
              trim = TRUE, 
              trim.cutoff = 0.015)

base_seq_ids = character(length(read_sets$consensus.sequences))
i = 1
for (seq_name in names(read_sets$consensus.sequences)){
  seq_id <- basename(seq_name)
  print(seq_id)
  new_seq_id <- paste(strsplit(seq_id, split_char)[[1]][[1]], seq_id)
  base_seq_ids[i] <- new_seq_id
  i <- i + 1
}

names(read_sets$consensus.sequences) <- base_seq_ids

read_sets$consensus.summaries
#read_sets$consensus.sequences
BrowseSeqs(read_sets$consensus.sequences)
 
fname = paste0(basename(input_folder),"_reads",".fasta")
setwd(args$output_folder)
write.dna(read_sets$consensus.sequences,
            file = fname,
           format = 'fasta',
           nbcol = -1,
           colsep = "",
           colw = 10000000)