#! /usr/bin/env python 
from Bio import Entrez
import argparse
import pprint
import sys


__author__ = "Andrew Ndhlovu"
__copyright__ = "Copyright 2018"
__license__ = "GPL"
__version__ = "3"
__maintainer__ = "Andrew Ndhlovu"


Entrez.email = "drewxdvst@outlook.com" 
query ="Phaeodactylum tricornutum"

def get_taxid(q_list):
    
    

    
    if q_list:
        for query in q_list:
            query = query.strip()
            record = esearch(query)
            if record:
                to_stdout(query, record)
            
    else:
        for query in sys.stdin:
            query = query.strip()
            record = esearch(query)
            if record:
                to_stdout(query, record)
            
            
def esearch(query):

    get  =  lambda query: Entrez.esearch(db="taxonomy", term =  query + "[Scientific Name]")
    
    handle = get(query)
    i  = 0
    while not handle:
        handle = get(query)
        i += 1
        if i == 5:
            sys.stderr.write(' '.join([query, "Fail\n"]))
            return
        
    return Entrez.read(handle)



def to_stdout(query, record):

    if len(record['IdList']) != 0:
        taxid = record['IdList'][0]
        print("{}\t{}".format(query, taxid))
        
    else:
        out = record['WarningList']['OutputMessage'][0] + "\n"
        sys.stderr.write("{}\t{}".format(query, out))
        
        
    

                       
        

if __name__ == '__main__':        
    parser = argparse.ArgumentParser(description="Retrieve a taxid for taxonomic name")
    group = parser.add_mutually_exclusive_group(required=False)
    group.add_argument('-q','--queries', action='store',nargs='+', type=str, required=False,
                        help='space seperated list of quoted queries ')
    group.add_argument('-l','--list-file', dest='list_file', action='store', required=False, type=str,
                        help='files with list of accessions ids')
    args = parser.parse_args()
    
    if args.queries:
        q_list = args.queries
    elif args.list_file:
        q_list = open(args.list_file.read().splitlines())
    else:
        q_list = None 
    get_taxid(q_list)
