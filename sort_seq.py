#! /usr/bin/env python
from Bio import SeqIO
import argparse
import pprint
import os

__author__ = "Andrew Ndhlovu"
__copyright__ = "Copyright 2018"
__license__ = "GPL"
__version__ = "3"
__maintainer__ = "Andrew Ndhlovu"
__email__ = "drewxdvst@outlook.com"


class Sort(object):

    def __init__(self):

        parser = argparse.ArgumentParser(description = """Sort sequences by id""")
        parser.add_argument('sequences',action='store', type=str)
        parser.add_argument('-i','--format', dest='ifmt', action='store', type=str, required=False,
                        default = "fasta",  help='format of input file')
        parser.add_argument('-o','--organism-name', dest = 'name_sort', action = 'store_true', help='sort by organism name')
        parser.add_argument('-f','--out-format', dest='ofmt', action='store', type=str, required=False,
                        help='format of output file')
        
        args  =  parser.parse_args()
        self.sequences = args.sequences
        self.name_sort = args.name_sort
        self.sorted = os.path.join(os.path.dirname(args.sequences),
                                    'sorted_'+os.path.basename(args.sequences))
        self.ifmt = args.ifmt
        self.ofmt = self.ifmt if not args.ofmt else args.ofmt  
        if self.name_sort and self.ifmt == 'fasta':
            print('Invalid argument, cannot sort using orginism name in fasta file, provide Genbank file')
            parser.print_help()

            
    def sort(self):

        sorted_seqs = []
        seqs =  SeqIO.parse(self.sequences, self.ifmt)
        ref_dict, seq_dict ={}, {}
        
        for n, rec in enumerate(seqs):
            if self.name_sort:
                ref_dict[n]=rec.annotations['organism']
            else:
                ref_dict[n]=rec.id
    
            seq_dict[n]=rec
            
        sorted_seqs = [ seq_dict[k] for k,v in sorted(iter(ref_dict.items()),\
                                                      key=lambda k_v: k_v[1])]
        with open(self.sorted, 'w') as fp:
            SeqIO.write(sorted_seqs, fp, self.ofmt)
            
if __name__ == '__main__':
    sort = Sort()
    sort.sort()
    
