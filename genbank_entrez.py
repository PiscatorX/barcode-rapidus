#! /usr/bin/env python 

from Bio import Entrez
from Bio import SeqIO
import argparse
import sys

__author__ = "Andrew Ndhlovu"
__copyright__ = "Copyright 2018"
__license__ = "GPL"
__version__ = "3"
__maintainer__ = "Andrew Ndhlovu"
__email__ = "drewxdvst@outlook.com"
Entrez.email = "drewxdvst@outlook.com"


def  GB_get(args):

    if args.list_file:
        seq_list = [ seq_id.split()[0] for seq_id in open(args.list_file).read().splitlines() ]
    if args.acc_ids:
        seq_list = args.acc_ids

    output_fname =  '.'.join([args.output_fname, args.out_format])
    f_obj =  open(output_fname,'w')
    
    for seq_id in seq_list:
        print seq_id
        handle = False
        while not handle:
            handle = Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id=seq_id)                            
            if handle:
                f_obj.write(handle.read())
                f_obj.flush()
                
    f_obj.close()


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description="Retrieve a list of genbank accession gis")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-l','--list-file', dest='list_file', action='store', required=False, type=str,
                        help='files with list of accessions ids')
    group.add_argument('-a','--acc-ids',dest='acc_ids', action='store',nargs='+', type=str, required=False,
                        help='space seperated list of  accession ids')
    parser.add_argument('-f','--out-format', dest='out_format', action='store', default='gb', required=False, type=str)
    parser.add_argument('-r','--retmode', dest='retmode', action='store', default='text', required=False, type=str)
    parser.add_argument('-o','--output-fname', dest='output_fname', action='store', default='GB_entrez', required=False, type=str)
    GB_get(parser.parse_args())
