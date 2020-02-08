#! /usr/bin/env python 

from xml.dom.minidom import parse
import urllib
from Bio import Entrez
import argparse
import time
import sys
import xml

__author__ = "Andrew Ndhlovu"
__copyright__ = "Copyright 2018"
__license__ = "GPL"
__version__ = "3"
__maintainer__ = "Andrew Ndhlovu"
__email__ = "drewxdvst@outlook.com"



def  GB_get(args):

    if args.list_file:
        seq_list = set([ seq_id.rstrip().split()[0] for seq_id in open(args.list_file).read().splitlines() ])
    if args.acc_ids:
        seq_list = args.acc_ids
    
    for acc_id in seq_list:
        query = "http://www.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=gene&id={}&rettype=text".format(acc_id)
        response = urllib.urlopen(query).read()
        print(response)
        #doc  = Entrez.parse(response, validate=False)
        #print(doc)
        #Idlist=map(getText, doc.getElementsByTagName('Id'))


        
def getText(node):
    assert node.firstChild.nodeType == node.TEXT_NODE
    return node.firstChild.data

        

        
        
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description="Retrieve a list of genbank accession gis")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-l','--list-file', dest='list_file', action='store', required=False, type=str,
                        help='files with list of accessions ids')
    group.add_argument('-a','--acc-ids',dest='acc_ids', action='store',nargs='+', type=str, required=False,
                        help='space seperated list of  accession ids')
    GB_get(parser.parse_args())




# idlist=map(_getText, doc.getElementsByTagName('Id'))
