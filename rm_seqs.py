#! /usr/bin/env python
from Bio import SeqIO
import argparse
import pprint
import  sys


def get_seqs(args):
        
    to_save = []
    rm_list = False
    
    if args.rm_list: 
       rm_list = open(args.rm_list).read().splitlines()
    if args.acc_list:
       rm_list = args.acc_list
                  
    with open(args.in_file) as fp:
       data = SeqIO.parse(fp, args.in_format)
       for seq in data:
           ref = seq.id
           if ref  in rm_list:
               print 'removed ',seq.id
               continue
           to_save.append(seq)
        
           
    if not args.out_format:
        args.out_format =  args.in_format
    
    if args.out_file:
        out_file  = args.out_file
    else:
        out_file = 'pruned_'+args.in_file.replace(args.in_format, args.out_format)
        
    if any(to_save):
        SeqIO.write(to_save, open(out_file,'w'), args.out_format)
    
            
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="""remove sequences from file""")
    parser.add_argument('in_file', action='store', type=str)
    parser.add_argument('-o','--out-file', dest='out_file', 
                        action='store', type=str)    
    parser.add_argument('-i','--in-format', dest='in_format', 
                        action='store', default="fasta", required=False, type=str)
    parser.add_argument('-f','--out-format', dest='out_format', 
                        action='store', required=False, type=str)
    group = parser.add_mutually_exclusive_group(required=True)    
    group.add_argument('-r','--remove-list', dest='rm_list', action='store', required=False, type=str)
    group.add_argument('-a','--acc-list', dest='acc_list', action='store', nargs='+', required=False, type=str)
    args = parser.parse_args()
    get_seqs(args)
