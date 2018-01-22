#!/usr/bin/python
from Bio import SeqIO
import pprint
import argparse
import sys
from collections import Counter, OrderedDict 

# seq_data = {}
# for line in open(sys.argv[1]).read().splitlines():
#     k, v = line.split("  ",1)
#     if k in 'ACGT':
#         continue
#     seq_data[k] = v.strip()
    
# pprint.pprint(seq_data)
# def getNs(seq_file):
#     with open(seq_file) as seq_fp:
#         records = SeqIO.parse(seq_fp, 'fasta')
#         print "="*50
#         print "{}\t{}\t{}\t{}".format("Len", "N_count", "N_perc", "seq_id")
#         print "="*50
#         for rec in records:
#             seq = rec.seq
#             seq_len, Ns = len(seq), seq.count('N')
#             try:
#                 value = round(float(Ns)/seq_len*100,2)
#             except:
#                 value = 0.00
#             print "{}\t{}\t{}%\t{}".format(seq_len, Ns, value, rec.id)

            
# getNs(sys.argv[1])

class CountAmbigousDNA(object):
    
    def __init__(self):

        
        parser = argparse.ArgumentParser(description="Counts numbers of ambiguous DNA in sequences",
                                         epilog="Requires Biopython")
        parser.add_argument('sequence_file', action='store')
        parser.add_argument('-f','--format', action='store', default='fasta', required=False,type=str)
        args = parser.parse_args()
        self.sequence_file = args.sequence_file
        self.format =  args.format
        self.ambiguous_DNA = ['B', 'D', 'H', 'K', 'M', 'N', 'S', 'R', 'U', 'W', 'V', 'Y', 'X']
        

    def  CountAmbiguousDNA(self):

         with open(self.sequence_file)  as fp:
             seq_data  = SeqIO.parse(fp, self.format)
             print '='*130
             print ''.join([ "{0:<10}".format(key) for key in self.ambiguous_DNA ])
             print '='*130
             for seq in seq_data:
                 n  = len(seq)
                 counter_seq  =  Counter(seq.seq)          
                 self.counter(n, counter_seq)
             print '='*130
             


             
    def counter(self, n, counter_seq):
        ambiguous_counts = []
        for code in self.ambiguous_DNA:
            count = counter_seq.get(code, 0)
            fraction  =  0 if count == 0 else count/float(n)
            ambiguous_counts.append("{0:}|{1:<8.1%}".format(count,fraction,code))
        print ''.join(ambiguous_counts)

            
if  __name__ == '__main__':
    sequence_data  = CountAmbigousDNA()
    sequence_data.CountAmbiguousDNA()
    print"""\n\n    Code Ambiguity  Complement Mnemonic
    A    A          T          adenine
    B    CGT        V          not_adenine
    C    C          G          cytosine
    D    AGT        H          not_cytosine
    G    G          C          guanine
    H    ACT        D          not_guanine
    K    GT         M          keto_base
    M    AC         K          amino_base
    N    ACGT       N          any_base
    R    AG         Y          purine_base
    S    CG         S          strong_bond
    T    T          A          thymine
    U    T          A          uracil
    V    ACG        B          not_thymine/uracil
    W    AT         W          weak_bond
    X    ACGT       X          unknown
    Y    CT         R          pyrimidine\n\n"""

    