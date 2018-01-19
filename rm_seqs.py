#! /usr/bin/env python
from Bio import SeqIO
import argparse
import pprint
import  sys


def get_seqs(args):
        
    to_save = []
    rm_list = False
    rm_to_save = []
    
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
               rm_to_save.append(seq)
               continue
           to_save.append(seq)
    
        
    if any(to_save):
        if not args.out_format:
            args.out_format =  args.in_format
        if args.out_file:
            out_file  = args.out_file
        else:
            out_file = 'pruned_'+args.in_file.replace(args.in_format, args.out_format)
            SeqIO.write(to_save, open(out_file,'w'), args.out_format)

    
    if args.rm_save and any(rm_to_save):
        if args.rm_out_file:
            out_file  = args.rm_out_file
        else:
            out_file = 'removed_'+args.in_file.replace(args.in_format, args.out_format)
        
        
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="""remove sequences from file using accession identifiers""")
    parser.add_argument('in_file', action='store', type=str)
    parser.add_argument('-o','--out-file', dest='out_file', 
                        action='store', type=str)    
    parser.add_argument('-i','--in-format', dest='in_format', 
                        action='store', default="fasta", required=False, type=str)
    parser.add_argument('-f','--out-format', dest='out_format', 
                        action='store', required=False, type=str)
    parser.add_argument('-s','--save-removed', dest='rm_save', 
                        action='store', required=False, type=str)
    parser.add_argument('-R','--rm-outfile', dest='rm_out_file', 
                        action='store', required=False, type=str)
    group = parser.add_mutually_exclusive_group(required=True)    
    group.add_argument('-r','--remove-list', dest='rm_list', action='store', required=False, type=str)
    group.add_argument('-a','--acc-list', dest='acc_list', action='store', nargs='+', required=False, type=str)
    args = parser.parse_args()
    if args.rm_out_file:
        args.rm_save = True
    get_seqs(args)
