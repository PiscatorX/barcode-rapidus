#!/usr/bin/env python 
from Bio import SeqIO
import argparse
import os

__author__ = "Andrew Ndhlovu"
__copyright__ = "Copyright 2019"
__license__ = "GPL"
__version__ = "3"
__maintainer__ = "Andrew Ndhlovu"
__email__ = "drewxdvst@outlook.com"



def extract_select():

    parser = argparse.ArgumentParser(description="Get select sequences from aligned/unaligned sequence. This utility does not validate input sequence format. Incase of empty output check input sequence format.", epilog="Requires Biopython")
    parser.add_argument('infile', action='store', type=str)
    parser.add_argument('-f','--format', dest='in_format', action='store', default='fasta', required=False, type=str)
    parser.add_argument('-o','--out-format', dest='out_format', action='store', required=False, type=str)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-l','--list-file', dest='list_file', action='store', required=False, type=str,
                        help='files with list of accessions ids')
    group.add_argument('-a','--acc-ids',dest='acc_ids', action='store',nargs='+', type=str, required=False,
                        help='space seperated list of  accession ids')
    parser.add_argument('-x','--extract-fname', dest='extract_fname', action='store', type=str, required=False,
                        help='output filename')
    
    args =  parser.parse_args()
    
    if  args.list_file:
        extract_list = [x.strip() for x in open(args.list_file).read().splitlines()]
    if  args.acc_ids:
        extract_list = args.acc_ids
        
    seqs = SeqIO.parse(open(args.infile), args.in_format)

    if not args.out_format:
        args.out_format = args.in_format
        
    if args.extract_fname:
        new_filename  = args.extract_fname
    else:
        path, fname = os.path.split(args.infile)
        filename, ext = os.path.splitext(fname)
        new_filename = 'extract_' + '.'.join([filename, args.out_format])

    SeqIO.write([ x for x in seqs if x.id in extract_list ],open(new_filename, 'w'), args.out_format)
    

if  __name__ == '__main__':
    extract_select()
