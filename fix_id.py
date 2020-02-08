#!/usr/bin/python


__author__ = "Andrew Ndhlovu"
__copyright__ = "Copyright 2018"
__license__ = "GPL"
__version__ = "3"
__maintainer__ = "Andrew Ndhlovu"
__email__ = "drewxdvst@outlook.com"



import argparse
from Bio import SeqIO
import sys
import os

def fix_id(fasta_file, spliton):

    fixed_fasta = os.path.join(os.path.dirname(fasta_file),
                                    'fixed_'+os.path.basename(fasta_file))
    fixed_fp = open(fixed_fasta,'w')
    with open(fasta_file) as fasta_object:
        fasta_seqs  = SeqIO.parse(fasta_object, 'fasta')
        new_seqs = []
        for i,seq  in enumerate(fasta_seqs, 1):
            split_id = seq.id.split(spliton,1)[0]
            #[0].split('-')[0]
            print(i,split_id)
            seq.id = split_id
            new_seqs.append(seq)
        SeqIO.write(new_seqs, fixed_fp, 'fasta')
    fixed_fp.close()
if  __name__ == '__main__':
    parser = argparse.ArgumentParser(description = """fix sequences id by splitting on provide char""")
    parser.add_argument('fasta_file',action='store', type=str, help="sequence in fasta format")
    parser.add_argument('-s','--spliton', default = '_', type=str, required=False,
                        help='split char for sequence id')
    args  =  parser.parse_args()
    fix_id(args.fasta_file, args.spliton)
