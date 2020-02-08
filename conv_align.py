#!/usr/bin/env python
from Bio import AlignIO
from Bio.Alphabet import generic_dna
import argparse
import os

__author__ = "Andrew Ndhlovu"
__copyright__ = "Copyright 2018"
__license__ = "GPL"
__version__ = "3"
__maintainer__ = "Andrew Ndhlovu"
__email__ = "drewxdvst@outlook.com"



def conv(args):


    if not args.outfile:
         outfile = ''.join([args.infile.replace(args.infmt, ''), args.ofmt])
    else:
        outfile  = args.outfile 
         
    with open(args.infile, "r") as input_handle, open(outfile, "w") as output_handle:
        AlignIO.convert(input_handle, args.infmt, output_handle, args.ofmt, generic_dna)    
    print('[{} :=> {}] {}'.format(args.infmt, args.ofmt, outfile))
    
if  __name__ == '__main__':
    parser = argparse.ArgumentParser(description="simple Biopython format convertor, converts gb to fasta by default")
    parser.add_argument('infile', action='store', type=str)
    parser.add_argument('-i','--in-format', dest='infmt', action='store', default='gb', type=str)
    parser.add_argument('-f','--out-format', dest='ofmt', action='store', default='fasta', required=False, type=str)    
    parser.add_argument('-o','--outfile', dest='outfile', action='store', type=str)
    args = parser.parse_args()
    conv(args)
