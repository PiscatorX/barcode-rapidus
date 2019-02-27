#!/usr/bin/env python
from collections import Counter
from  Bio import SeqIO
import argparse
import sys



def parse_stats(f_name, fmt):
    msa = SeqIO.parse(f_name, fmt)
    i = 0
    ambiguous_DNA = ['B', 'D', 'H', 'K', 'M', 'N', 'S', 'R', 'U', 'W', 'V', 'Y', 'X']
    for seq in msa:
        i += 1
        n = len(seq)
        ambigs  = sum([ v for k,v in Counter(seq.seq).items() if k in ambiguous_DNA ])
        perc_ambigs = 0 if n == 0 else ambigs/float(n)
        if fmt=='gb':  
            print "{}\t{}\t{}({:.2%})\t{}\t{}".format(i, n, ambigs, perc_ambigs, seq.annotations['organism'], seq.id)
        else:
            print "{}\t{}\t{}({:.2%})\t{}".format(i, n, ambigs, perc_ambigs, seq.id)



if __name__ ==  '__main__':
    parser = argparse.ArgumentParser(description = """get sequence stats""")
    parser.add_argument('sequences',action='store', type=str)
    parser.add_argument('-f','--format', dest='format', 
                        action='store', default="fasta", required=False, type=str)
    args  =  parser.parse_args()
    
    parse_stats(args.sequences, args.format)
