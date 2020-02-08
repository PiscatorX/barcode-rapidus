#!/usr/bin/env python
import argparse
import collections
from Bio import SeqIO

def dedup(seq_fobj, informat, seq_out, outformat):
    
    fasta_records = SeqIO.parse(seq_fobj, informat)
    id_counts =  collections.defaultdict(int)
    fid_counts =  collections.defaultdict(int)
    duplicates = collections.defaultdict(list)
    dedup_i = 0
    
    for i,seq  in enumerate(fasta_records):
       id_counts[seq.id]+=1
       seq_id = seq.id+"_"+str(id_counts[seq.id]) if id_counts[seq.id] != 1 else seq.id
       if id_counts[seq.id] != 1:
           duplicates[seq.id].append(seq_id)
           print("{}\tfasta duplicate: {}\t{}".format(dedup_i, seq.id, seq_id))
           dedup_i += 1
           seq.id = seq_id
       SeqIO.write(seq, seq_out, outformat)
       seq_out.flush()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Dedupe fasta""")
    parser.add_argument('sequences', type=argparse.FileType('r') )
    parser.add_argument('-i','--informat', default="fasta")
    parser.add_argument('-o','--outfile', type=argparse.FileType('w'), default="deduped.out")
    parser.add_argument('-f','--format', default="fasta")
    args, unknown = parser.parse_known_args()
    dedup(args.sequences, args.informat, args.outfile, args.format)
    
