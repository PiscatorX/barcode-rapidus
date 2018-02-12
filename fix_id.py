#!/usr/bin/python

from Bio import SeqIO
import sys
import os

def fix_id(fasta_file):

    fixed_fasta = os.path.join(os.path.dirname(fasta_file),
                                    'fixed_'+os.path.basename(fasta_file))
    fixed_fp = open(fixed_fasta,'w')
    with open(fasta_file) as fasta_object:
        fasta_seqs  = SeqIO.parse(fasta_object, 'fasta')
        new_seqs = []
        for i,seq  in enumerate(fasta_seqs, 1):
            split_id = seq.id.split('_',1)[0]
            #[0].split('-')[0]
            print i,split_id
            seq.id = split_id
            new_seqs.append(seq)
        SeqIO.write(new_seqs, fixed_fp, 'fasta')
    fixed_fp.close()

fix_id(sys.argv[1])
