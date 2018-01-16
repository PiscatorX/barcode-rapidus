#!/usr/bin/python
import sys
import time
import argparse
from Bio import Entrez

__author__ = "Andrew Ndhlovu"
__copyright__ = "Copyright 2018"
__license__ = "GPL"
__version__ = "3"
__maintainer__ = "Andrew Ndhlovu"
__email__ = "drewxdvst@outlook.com"


Entrez.email = "drewxdvst@outlook.com" 


class  GenBankGet(object):

    def  __init__(self):

        parser = argparse.ArgumentParser(description="""Retrieve sequences from GenBank using query e.g "(Prorocentrum[Primary Organism] AND 28S[Text Word])" """,
                                         epilog="Requires Biopython")
        group = parser.add_mutually_exclusive_group(required=True)
        group.add_argument('-q','--query', dest='query', action='store', required=False, type=str)
        group.add_argument('-l','--query-file', dest='query_list', action='store', required=False,type=str,
                            help="list of genbank queries on each line of the file")
        parser.add_argument('-m','--retmax', dest='retmax', action='store', required=False, default=15,type=int)
        parser.add_argument('-b','--batchsize', dest='batchsize', action='store', required=False, default=5,type=int)

        
        self.args = parser.parse_args()
        if  self.args.query_list:
            self.query = open(self.args.query_list).read().splitlines()
        if  self.args.query:
            self.query = [self.args.query]

    def query_genbank(self):
        map(self.gb_search, self.query)

        
    def gb_search(self, query):

         with open('query_genbank_history','a') as fp:
             fp.writelines([query])
             fp.write('\n')
             fp.flush()
             
         print query
         fail = 0
         handle = 0
         while not handle:
             try:
                 handle = Entrez.esearch(db="Nucleotide", term=query, usehistory="y") 
             except Exception as e:
                 fail += 1
                 if fail == 15:
                     print "Failed to connect to NCBI entrez server\nExiting...\n"
                     sys.exit(1)
                 print e
                 time.sleep(5)
                 print "Reconnecting to NCBI entrez..."
                 continue

         record = Entrez.read(handle)
         n_rec =  int(record['Count'])
         if self.args.retmax:
             n_rec = self.args.retmax
         self.args.batchsize =  n_rec if self.args.batchsize >  n_rec else n_rec  
         
         with open('genbank_'+time.strftime('%A-%H-%M-%S')+'.gb','w') as fp:
             print
             for retstart  in range(0, n_rec, self.args.batchsize):
                 handle = False
                 fail = 0
                 while not handle:
                     try:
                         handle = Entrez.efetch(db="Nucleotide", rettype="gb",retmode="text", 
                          retstart=retstart,retmax=self.args.retmax, webenv=record['WebEnv'], 
                          query_key=record["QueryKey"])
                         
                     except Exception as e:
                         fail += 1
                         if fail == 15:
                             print "Failed to connect to NCBI entrez server\nExiting...\n"
                             sys.exit(1)
                         print e
                         time.sleep(5)
                         print "Reconnecting to NCBI entrez..."
                         continue

                 results =  handle.read()
                 fp.write(results)
                 count = retstart+self.args.batchsize
                 if count > n_rec:
                     count = n_rec
                     
                 args='\r{0}/{1} ({2:.2%})\r'.format(count, n_rec,  (count)/float(n_rec))
                 sys.stdout.write(args)
                 sys.stdout.flush()
                 time.sleep(2)
                 
if __name__ == '__main__':
     gb = GenBankGet()
     gb.query_genbank()
