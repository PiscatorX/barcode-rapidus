#! /usr/bin/env python
from Bio import SeqIO
import argparse

__author__ = "Andrew Ndhlovu"
__copyright__ = "Copyright 2018"
__license__ = "GPL"
__version__ = "3"
__maintainer__ = "Andrew Ndhlovu"
__email__ = "drewxdvst@outlook.com"



def cli_init():

    parser = argparse.ArgumentParser(description="remove duplicate sequences and save list of seq ids", epilog="Requires Biopython")

    parser.add_argument('infile', action='store', type=str)
    parser.add_argument('-f','--format', dest='in_format', action='store', default='gb', required=False, type=str)
    parser.add_argument('-o','--out_format', dest='out_format', action='store', default='gb', required=False, type=str)
    return parser.parse_args()


def rm_dup(infile, in_format, out_format):

    seqs = SeqIO.parse(open(infile), in_format)
    duduped_seqs = []
    seq_ids = []
    for rec  in seqs:
        id_tmp  = rec.id
        if id_tmp in seq_ids:
            print "removed {}.".format(id_tmp)
            continue
        
        seq_ids.append(id_tmp)
        duduped_seqs.append(rec)
    seq_ids.sort()
    out_f = 'dedup_'+infile.replace(in_format, out_format)
    
    SeqIO.write(duduped_seqs, open(out_f,'w'), out_format)
    
    with open(infile+'.ids','w') as fp:
        fp.writelines(map(lambda x: x+'\n', seq_ids))

if  __name__ == '__main__':
    args = cli_init()
    rm_dup(args.infile, args.in_format, args.out_format)
