#! /usr/bin/env Rscript

setwd("~/Documents/barcode-rapidus")
require(ape)
require(reshape2)
require(phangorn)
require(stringi)
require(stringr)
require(argparse)
require(ggplot2)
require(here)
#install.packages('grid')
#require('grid')
#install.packages('gridExtra')
#require(gridExtra)
# Bioconductor packages
require(DECIPHER)
require(Biostrings)
require(sangerseqR)
#github("roblanf/sangeranalyseR")
require(sangeranalyseR)


parser <- ArgumentParser()
parser$add_argument("-i", "--input-folder",nargs='1',
                    help="input folder where abi files are located.")

args <- parser$parse_args()


# make.readsets(input.folder, forward.suffix, reverse.suffix, trim = TRUE,
#               trim.cutoff = 1e-04, processors = NULL, min.length = 1,
#               max.secondary.peaks = NULL, secondary.peak.ratio = 0.33)
# 
# wd =  cwd
# sf = summarise.abi.folder(wd)
# sf$summaries
# require(gplot2)
# ggplot(sf$summaries, aes(x = folder.name, y = raw.mean.quality)) + geom_boxplot()
# ggsave()
# 
# ggplot(sf$summaries, aes(x = folder.name, y = trimmed.mean.quality)) + geom_boxplot()
# ggsave()
# 
# # In this case, we also add a add horizontal lines at a cutoff of 2
# ggplot(sf$summaries, aes(x = folder.name, y = trimmed.secondary.peaks)) + geom_boxplot() + geom_hline(yintercept = 2, linetype = 3, colour = 'red')
# 
# ggplot(sf$summaries, aes(x = trimmed.mean.quality, y = trimmed.secondary.peaks)) + geom_point()



