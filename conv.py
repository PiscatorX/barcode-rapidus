#!/usr/bin/python
from Bio import AlignIO
from Bio.Alphabet import generic_dna
import argparse
import os




def conv(infile, infmt, outfile, ofmt):
    
    if not outfile:
         outfile = '.'.join([os.path.splitext(infile)[0], ofmt])
         print 'Converting to {}'.format(outfile)

    with open(infile, "rU") as input_handle, open(outfile, "w") as output_handle:  
        alignments = AlignIO.parse(input_handle, infmt, generic_dna)
        AlignIO.write(alignments, output_handle, ofmt)
    
    
if  __name__ == '__main__':
    parser = argparse.ArgumentParser(description="format convertor", epilog='used the Biopython modules')
    parser.add_argument('infile', action='store', required=True, type=str)
    parser.add_argument('-i','--in-format', dest='infmt', action='store', default='fasta', type=str)
    parser.add_argument('-f','--out-format', dest='ofmt', action='store', required=True, type=str)    
    parser.add_argument('-o','--outfile', dest='outfile', action='store', type=str)
    return parser.parse_args()

args = get_args()
conv(args.infile, args.infmt, args.outfile, args.ofmt)
