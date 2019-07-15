#! /usr/bin/env python
from Bio import SeqIO
import argparse

__author__ = "Andrew Ndhlovu"
__copyright__ = "Copyright 2018"
__license__ = "GPL"
__version__ = "3"
__maintainer__ = "Andrew Ndhlovu"
__email__ = "drewxdvst@outlook.com"


def show_ids():
    
    parser = argparse.ArgumentParser(description="save list of seq ids", epilog="Requires Biopython")
    parser.add_argument('infile', action='store', type=str)
    parser.add_argument('-f','--format', dest='in_format', action='store', default='gb', required=False, type=str)
    args = parser.parse_args()

    with open(args.infile+'.ids','w') as fp:
        seqs = SeqIO.parse(open(args.infile), args.in_format)
        for rec  in seqs:
            print('# '+rec.id, file=fp)
            print(rec.id)
    
if  __name__ == '__main__':
    show_ids()
