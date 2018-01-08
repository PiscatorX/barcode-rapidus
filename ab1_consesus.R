#! /usr/bin/env Rscript

require(ape)
require(reshape2)
require(phangorn)
require(stringi)
require(stringr)
#(devtools)
require(argparse)
# Bioconductor packages
require(DECIPHER)
require(Biostrings)
require(sangerseqR)
#github("roblanf/sangeranalyseR")
require(sangeranalyseR)

parser <- ArgumentParser()
parser$add_argument("abi",nargs='+', help="list of ab1 files")
args <- parser$parse_args()

print(args$abi)
# abi_file1="/home/drewx/Documents/sanpccX/DevOps/AlignDevX/raw/FSA5.1_rbcL-1385R_A02_02.ab1"
# abi_file2="/home/drewx/Documents/sanpccX/DevOps/AlignDevX/raw/FSA5.1_rbcL1385R_C12_09.ab1"
#             seq_filepath = args_abi
# seq_abif = read.abif(seq_filepath)
# seq_sanger = sangerseq(seq_abif)
# 
# trims = trim.mott(seq.abif, cutoff=0.015)
# seq_untrimmed = seq.abif@data$PBAS.2
# seq_trimmed = substring(seq.untrimmed, trims$start, trims$finish)
# print(seq_trimmed)
# print(paste("Untrimmed Sequence length:", nchar(seq.abif@data$PBAS.2)))
# print(paste("Trimmed Sequence length:", nchar(seq_trimmed)))            




