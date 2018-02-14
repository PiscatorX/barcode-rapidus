#! /usr/bin/env  python
from Bio import  SeqIO
import  argparse

__author__ = "Andrew Ndhlovu"
__copyright__ = "Copyright 2018"
__license__ = "GPL"
__version__ = "3"
__maintainer__ = "Andrew Ndhlovu"
__email__ = "drewxdvst@outlook.com"




class  Features(object):

    def __init__(self):

        parser = argparse.ArgumentParser(description="""retrive features from a gb file. Sequence annotation information is preseved. Once edited all feature positional information must be assumed to incorrect for those sequences where features have been extracted.""")    
        parser.add_argument('genbank_file', action='store', type=str)
        parser.add_argument('-k','--key', dest='key', action='store', required=False, default='CDS', type=str,
                            help="This is a textual description of the type of feature e.g. 'CDS' or 'gene'")
        parser.add_argument('-q','--qualifier', dest='qualifier', action='store', required=False, default='gene', type=str)
        parser.add_argument('-n','--q_name', dest='q_name', action='store', required=False, default='rbcL', type=str)
        args = parser.parse_args()
        self.genbank_file = args.genbank_file
        self.key = args.key
        self.qualifier = args.qualifier
        self.q_name = args.q_name

        
    def parse_gb(self):
        
        seq_records = SeqIO.parse(self.genbank_file, "gb")
        self.new_seq_records = []
        count =  0
        for record in seq_records:
            
            m = len(record)
            extract_feature = self.get_feature(record)
            annotations = record.annotations
            if not extract_feature:
                 self.new_seq_records.append(record)
            else:
                 n = len(extract_feature)
                 if (m == n):
                     self.new_seq_records.append(record)
                 else:
                     #print extract_feature
                     new_record = extract_feature.extract(record)
                     new_record.annotations = annotations
                     self.new_seq_records.append(new_record)
                     print '>>> {}\t{}\t{}\t{}'.format(record.id, record.annotations['organism'], m, n)
                     count += 1
        print '\n{} features extracted!'.format(count)
    def get_feature(self, record):
        
        for feature in record.features:
            try:
                if feature.type == self.key and feature.qualifiers[self.qualifier] == [self.q_name]:
                   target = feature.location
                   if target.strand == 1:
                       return target
            except KeyError:
                return False
            
        return False

    

    def save_fx(self):

        
        with open('fx_'+self.genbank_file,'w') as fp:
            SeqIO.write(self.new_seq_records, fp, 'gb')
        

if __name__ == '__main__':
        features = Features()
        features.parse_gb()
        features.save_fx()
