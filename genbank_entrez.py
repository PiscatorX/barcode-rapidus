#! /usr/bin/env python 

from Bio import Entrez
from Bio import SeqIO
import argparse
import time
import sys

__author__ = "Andrew Ndhlovu"
__copyright__ = "Copyright 2018"
__license__ = "GPL"
__version__ = "3"
__maintainer__ = "Andrew Ndhlovu"
__email__ = "drewxdvst@outlook.com"
Entrez.email = "drewxdvst@outlook.com"


def  GB_get(args, hist_file = "genebank.entrez.hist"):

    if args.list_file:
        seq_list = set([ seq_id.rstrip().split()[0] for seq_id in open(args.list_file).read().splitlines() ])
    if args.acc_ids:
        seq_list = args.acc_ids
    output_fname =  '.'.join([args.output_fname, args.out_format])
    f_obj =  open(output_fname,'w')

    if args.cont:
        if os.path.exists(hist_file):
            hist_list = [ seq_id.rstrip().split()[0] for seq_id in open(hist_file).read().splitlines() ]
            seq_list = set(seq_list) -  set(hist_list)
        else:
            print("Nothing to resume, '{}' file not found in currend directory".format(hist_file))
    
    hist_fobj = open(hist_file, "a")
    for seq_id in seq_list:
        print(repr(seq_id))
        fail = 0
        handle = False
        while not handle:
            try:
                handle = Entrez.efetch(db = args.db,
                                       id = seq_id,
                                       rettype = args.out_format,
                                       retmode ="text")
                if handle:
                    f_obj.write(handle.read())
                    f_obj.flush()
                    hist_fobj.write(seq_id+"\n")
                    hist_fobj.flush()
                    
            except Exception as e:
                print(e)
                fail += 1
                if 'Bad Request' in str(e):
                    break
                time.sleep(1)                                                                                                                                                       
                if (fail == 10):
                    sys.exit(1)
                print("Reconnecting to NCBI entrez. Re-try: ",fail)
                
    f_obj.close()
    hist_fobj.close()
     
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description="Retrieve a list of genbank accession gis")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-l','--list-file', dest='list_file', action='store', required=False, type=str,
                        help='files with list of accessions ids')
    group.add_argument('-a','--acc-ids',dest='acc_ids', action='store',nargs='+', type=str, required=False,
                        help='space seperated list of  accession ids')
    parser.add_argument('-f','--out-format', dest='out_format', action='store', default='gb', required=False, type=str)
    parser.add_argument('-r','--retmode', dest='retmode', action='store', default='text', required=False, type=str)
    parser.add_argument('-d','--db', default='nucleotide', choices = ['nucleotide', 'protein', 'nuccore' ], required=False, type=str)
    parser.add_argument('-o','--output-fname', dest='output_fname', action='store', default='GB_entrez', required=False, type=str)
    parser.add_argument('-c','--cont', help = "stores sucessfully downloaded ids in \"genebank.entrez.hist\" file ",  action='store_true')
    GB_get(parser.parse_args())
