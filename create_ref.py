#!/usr/bin/python

from Bio import  SeqIO
import argparse
import sys

def  get_ref(args):
     
     ref_obj = open(args.infile + '.ref','w')
     seq_obj = open('final_' + args.infile.replace(args.in_format,'') + args.out_format, 'w')
     with open(args.infile) as fp:
         seq_data = SeqIO.parse(fp, args.in_format)
         for i,record in enumerate(seq_data, 1):
             if  args.in_format == 'gb':
                  organism = record.annotations['organism']
             else:
                  organism = ' '.join(record.description.split()[1:3])
             refs ="{}\t{}".format(record.id, organism)     
             print >>ref_obj,refs
             print refs
             record.description = ''
             SeqIO.write(record, seq_obj , args.out_format)
     map(lambda f : f.close, [ref_obj,  seq_obj])

     
if  __name__ ==  '__main__':
    parser = argparse.ArgumentParser(description="creates a reference file of seq_id and organism names which", epilog="Requires Biopython")
    parser.add_argument('infile', action='store', type=str)
    parser.add_argument('-f','--format', dest='in_format', action='store', default='gb', required=False, type=str)
    parser.add_argument('-o','--out_format', dest='out_format', action='store', default='fasta', required=False, type=str)
    get_ref(parser.parse_args())

