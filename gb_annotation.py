#! /usr/bin/env python
from Bio import SeqIO


__author__ = "Andrew Ndhlovu"
__copyright__ = "Copyright 2018"
__license__ = "GPL"
__version__ = "3"
__maintainer__ = "Andrew Ndhlovu"
__email__ = "drewxdvst@outlook.com"



def data(args):

    anno = args.annotation 
    with open(args.in_file) as fp:
         data  = SeqIO.parse(fp, 'gb')
         for  count, seq  in enemerate(data,1):
                print count,seq.id,seq.annotations[anno]
    
if __name__  == '__main__':

    parser = argparse.ArgumentParser(description="""view and count annotations""")
    
    parser.add_argument('in_file', action='store', type=str)
    parser.add_argument('-a','--anno', dest='annotation', action='store', required=False, type=str)
    args = parser.args()
    get_seqs(args):
