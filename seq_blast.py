#! /usr/bin/python

from Bio.Blast import NCBIWWW
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




def get_sequences(query_sequences, fmt, output_dir, hitlist_size, max_queries):

    with open(query_sequences) as query_fp:
        seq_data = SeqIO.parse(query_fp, fmt)
        i = 0
        for seq in seq_data:
            if not seq.seq:
                continue
            result_id = '.'.join(['_' + str(i),'blast'])
            fname = ''.join([seq.id.replace('.ab1', ''), result_id])
            dest_dir = output_dir
            try:
                os.mkdir(dest_dir)
            except OSError:
                pass
            save_fname = os.path.join(dest_dir, fname)
            if  os.path.exists(save_fname):
                i += 1
                print i,save_fname,'File exists...'
                continue
            results_handle = None
            fail = 0
            while not results_handle:
                try:
                    results_handle = NCBIWWW.qblast("blastn", "nt", seq.format(fmt), hitlist_size=hitlist_size)
                except Exception as e:
                    fail += 1
                    if fail == 15:
                       print "Failed to connect to NCBIWWW server\nExiting...\n"
                       sys.exit(1)
                    print e
                    time.sleep(5)
                    print "Reconnecting to NCBIWWW server..."
                    continue
            with open(save_fname,'w') as save:
                 save.write(results_handle.read())
            i += 1
            print i,save_fname
            if i ==  max_queries:
                print "Max queries reached! Exiting..."
                sys.exit(0)
                
            time.sleep(5)
            

if __name__ ==  '__main__':
    
    parser = argparse.ArgumentParser(description="blast sequences against NCBIWWW online server")
    parser.add_argument('query_sequences', action='store', type=str,
                        help='sequenced or unsequenced query sequence')
    parser.add_argument('-f','--format', dest='fmt', action='store', type=str, required=False,
                        default = "fasta",  help='format of input file')
    parser.add_argument('-o','--ouput-dir', dest='output_dir', action='store', type=str, required=False,
                        default = "blast_results",  help='output directory')
    parser.add_argument('-s','--hitlist-size', dest='hitlist_size', action='store', type=int, required=False,
                        default =1,  help='output directory')
    parser.add_argument('-m','--max-queries', dest='max_queries', action='store', type=int, required=False,
                        default =100)
    args = parser.parse_args()
    get_sequences(args.query_sequences, args.fmt, args.output_dir, args.hitlist_size, args.max_queries)
                        
