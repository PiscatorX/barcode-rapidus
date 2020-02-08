#!/usr/bin/env python

from Bio import  SeqIO
import argparse
import sys
import os

__author__ = "Andrew Ndhlovu"
__copyright__ = "Copyright 2018"
__license__ = "GPL"
__version__ = "3"
__maintainer__ = "Andrew Ndhlovu"
__email__ = "drewxdvst@outlook.com"



def  get_ref(args):
     
     ref_obj = open(args.infile + '.ref','w')
     path, fname = os.path.split(args.infile)
     fname =  ''.join(['final_',os.path.splitext(fname)[0],'.', args.out_format])
     fname = os.path.join(path, fname)
     seq_obj = open(fname, 'w')
     with open(args.infile) as fp:
         seq_data = SeqIO.parse(fp, args.in_format)
         for i,record in enumerate(seq_data, 1):
             if args.uniprot:
                  organism = record.description.split("OS=")[-1].split("OX=")[0].strip()
                  record.id = record.id.split("|")[1]
             else: 
                  if  args.in_format == 'gb':
                      organism = record.annotations['organism']
                  else:
                      organism = ' '.join(record.description.split()[1:3])

             refs ="{}\t{}".format(record.id, organism)     
             print(refs, file=ref_obj)
             print(refs)
             record.description = ''
             SeqIO.write(record, seq_obj , args.out_format)
     list(map(lambda f : f.close, [ref_obj,  seq_obj]))

     
if  __name__ ==  '__main__':
    parser = argparse.ArgumentParser(description="creates a reference file of seq_id and organism names which", epilog="Requires Biopython")
    parser.add_argument('infile', action='store', type=str)
    parser.add_argument('-f','--format', dest='in_format', action='store', default='gb', required=False, type=str)
    parser.add_argument('-o','--out_format', dest='out_format', action='store', default='fasta', required=False, type=str)
    parser.add_argument('-u','--uniprot', action='store_true', default=False)
    get_ref(parser.parse_args())

