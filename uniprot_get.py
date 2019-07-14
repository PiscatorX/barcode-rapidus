#!/usr/bin/env python

import urllib,urllib2
from Bio import Entrez
from Bio import SeqIO
import argparse
import time
import sys
import os

__author__ = "Andrew Ndhlovu"
__copyright__ = "Copyright 2018"
__license__ = "GPL"
__version__ = "3"
__maintainer__ = "Andrew Ndhlovu"
__email__ = "drewxdvst@outlook.com"
contact = "drewxdvst@outlook.com"
url = "https://www.uniprot.org/uploadlists/"




def  uniprot_get(args, hist_file = "uniprot.hist"):

    if args.list_file:
        seq_list = set([ seq_id for seq_id in open(args.list_file).read().split() ])
    if args.acc_ids:
        seq_list = args.acc_ids
        
    output_fname =  '.'.join([args.outfname, args.format])
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
        handle  = False
        while not handle:
            try:
                params = {"from":"ACC+ID",
                          "to":"ACC",
                          "format":"fasta",
                          "query": seq_id }

                data = urllib.urlencode(params)
                request = urllib2.Request(url, data)
                request.add_header("User-Agent","Python %s" % contact)
                response = urllib2.urlopen(request)
                handle  = response.read()
                if not handle:
                    fail += 1
                    print("Reconnecting to uniprot. Re-try: ",fail)
                    if (fail == 5):
                        break
                    continue
                if handle:
                    f_obj.write(handle)
                    f_obj.flush()
                    hist_fobj.write(seq_id+"\n")
                    hist_fobj.flush()
                
                    
            except Exception as e:
                fail += 1
                print(e)                                                                                                                                
                if (fail == 10):
                    sys.exit(1)
                print("Reconnecting to uniprot. Re-try: ",fail)
                
            time.sleep(1) 
           
                            
    f_obj.close()
    hist_fobj.close()
     
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description="Retrieve a list of genbank accession gis")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-l','--list-file', dest='list_file', action='store', required=False, type=str,
                        help='files with list of accessions ids')
    group.add_argument('-a','--acc-ids',dest='acc_ids', action='store',nargs='+', type=str, required=False,
                        help='space seperated list of  accession ids')
    parser.add_argument('-f','--format', action='store', default='fasta', choices = ["fasta", "text"], required=False, type=str)
    parser.add_argument('-o','--outfname', action='store', default='uniprot_KB.files', required=False, type=str)
    parser.add_argument('-c','--cont', help = "stores sucessfully downloaded ids in \"uniprot.hist\" file ",  action='store_true')
    uniprot_get(parser.parse_args())
