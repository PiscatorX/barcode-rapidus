#! /usr/bin/env python
from Bio import SeqIO
from collections import defaultdict
import argparse
import pprint
import  sys

__author__ = "Andrew Ndhlovu"
__copyright__ = "Copyright 2018"
__license__ = "GPL"
__version__ = "3"
__maintainer__ = "Andrew Ndhlovu"
__email__ = "drewxdvst@outlook.com"


def get_seqs(args):
                  
    with open(args.in_file) as fp:
       data = SeqIO.parse(fp, args.in_format)
       seq_data =  defaultdict(list)   
       for seq in data:
           seq_data[seq.id].append(seq)
           
    for seq_id, seq_records in  list(seq_data.items()):
        n =  len(seq_records)
        if n <= 1:
            continue
        print("{}\t{}".format(seq_id, n))
        if args.save:
            seq_fname =  '.'.join([seq_id, args.out_format])
            with open(seq_fname, 'w') as fp:
                SeqIO.write(seq_data[seq_id], open(seq_fname,'w'), args.out_format)
                  
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="""identifies repeated sequences based on sequence identifiers""")
    parser.add_argument('in_file', action='store', type=str)
    parser.add_argument('-i','--in-format', dest='in_format', action='store', default="fasta", required=False, type=str)
    parser.add_argument('-f','--out-format', dest='out_format', action='store')
    parser.add_argument('-s','--save', dest='save', help="Save repeated sequences based on sequence identifier", action='store_true')
    args = parser.parse_args()
    if not args.out_format:
         args.out_format = args.in_format
    get_seqs(args)
