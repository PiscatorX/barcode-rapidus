#! /usr/bin/env python
from Bio import SeqIO
import argparse
import sys

__author__ = "Andrew Ndhlovu"
__copyright__ = "Copyright 2018"
__license__ = "GPL"
__version__ = "3"
__maintainer__ = "Andrew Ndhlovu"
__email__ = "drewxdvst@outlook.com"



def get_anno(args):

    anno = args.annotation 
    with open(args.genbank_file) as fp:
         data  = SeqIO.parse(fp, 'gb')
         for  count, seq  in enumerate(data,1):
                if args.description:
                    print(seq.id, seq.description)
                    continue
                sys.stderr.write(seq.id+'\t')
                print(seq.annotations[anno])
                
if __name__  == '__main__':

    parser = argparse.ArgumentParser(description="""convert and remove sequences""")
    
    parser.add_argument('genbank_file', action='store', type=str)
    group = parser.add_mutually_exclusive_group(required=False)
    group.add_argument('-d','--description', dest='description', action='store_true', required=False)
    group.add_argument('-a','--anno', dest='annotation', action='store', required=False, default='taxonomy', type=str)
    args = parser.parse_args()
    get_anno(args)
