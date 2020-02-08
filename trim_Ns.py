#! /usr/bin/env python
from Bio import SeqIO
import argparse


__author__ = "Andrew Ndhlovu"
__copyright__ = "Copyright 2018"
__license__ = "GPL"
__version__ = "3"
__maintainer__ = "Andrew Ndhlovu"
__email__ = "drewxdvst@outlook.com"



class Trim(object):

    def __init__(self):

        parser = argparse.ArgumentParser(description = """Trimm trailing and leading Ns from fastA file""")
        parser.add_argument('sequences',action='store', type=str)
        args  =  parser.parse_args()
        self.sequences = args.sequences
        self.trimmed = 'trimmed_'+args.sequences

        
    def trim(self):

        trimmed_seqs = []
        seqs =  SeqIO.parse(self.sequences, 'fasta')
        for rec in seqs:
            trimmed_seq = self.rm_N(rec.seq.upper())
            rec.seq = trimmed_seq
            if  not any(rec.seq):
                print("{} removed, no data!".format(rec.id))
                continue
            
            trimmed_seqs.append(rec)
            
        with open(self.trimmed, 'w') as fp:
            SeqIO.write(trimmed_seqs, fp, 'fasta')
        

    def rm_N(self, seq):

        lead = trail =  0
        if seq.startswith('N'):
             for lead, base in enumerate(seq):
                 if base != 'N':
                    break       
        if seq.endswith('N'):
            for trail, base in enumerate(seq[::-1]):
                if base != 'N':
                    break
        else:
            return seq[lead:]
                    
        return seq[lead:-trail] 
                

if __name__ == '__main__':
    trim = Trim()
    trim.trim()
