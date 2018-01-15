#!/usr/bin/python

from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO
import argparse
import sqlite3
import glob
import csv
import os
__author__ = "Andrew Ndhlovu"
__copyright__ = "Copyright 2018"
__credits__ = "Rob Knight"
__license__ = "GPL"
__version__ = "3"
__maintainer__ = "Andrew Ndhlovu"
__email__ = "drewxdvst@outlook.com"




class ParseBlast(object):

    def __init__(self):

        parser = argparse.ArgumentParser(description="parse blast results and optionally create a database of blast results")
        parser.add_argument('-r', '--results-dir', dest='results_dir', action='store', type=str, default="blast_results",
                            help="directory where blast files are located")
        parser.add_argument('-f','--fname', dest='fname_patt', action='store', required=False, default= "*.blast", type=str,
                            help="pattern of filenames/file extrensions to be matched")
        parser.add_argument('-l','--hit-limit', dest='hit_limit', action='store', required=False, default= 5, type=int,
                            help="shortlist of accession numbers per hit to be saved to file")
        parser.add_argument('-d','--db-name', dest='db_name',  action='store', required=False, type=str)

        args  =  parser.parse_args()
        file_patt =  args.fname_patt
        self.results_dir = args.results_dir
        self.blast_files = glob.iglob(os.path.join(os.path.join(self.results_dir, file_patt)))
        self.hit_limit = args.hit_limit
        self.drop_tags = ['sbjct', 'frame', 'query', 'strand', 'num_alignments', 'match']
        self.DB = args.db_name
        if self.DB:
             self.DB_conn  = sqlite3.connect(self.DB)
             self.DB_c = self.DB_conn.cursor()
        
    def get_data(self):
        
       fname_prefix = os.path.basename(self.results_dir)
       self.refs = {}
       
       with open(fname_prefix+'_blast.data','w') as data_file:
           for blast_result in self.blast_files:
               self.analyse_result(blast_result, data_file)
               if self.DB:
                   self.DB_conn.commit()
               
       with open(fname_prefix+'_blast.ids','w') as gi_file:
           for gi, acc in self.refs.items():
               print >>gi_file, "{}\t{}".format(gi, acc)

               
    def analyse_result(self, result_file, data_file):

       
       with open(result_file) as file_obj:
           header = '# '+result_file
           print header
           print >>data_file,header
           data = NCBIXML.parse(file_obj)
           records = data.next()
           query = records.query
           hsp_results = {'seq_id' : query }
           
           for count, hit in enumerate(records.alignments, 1):
               #print vars(hit)
               hsp = vars(hit.hsps[0])
               map(hsp.pop, self.drop_tags)
               map(hsp_results.update, [hsp, vars(hit)])
               hit_def = hsp_results['hit_def']
               hit_id = hsp_results['hit_id']
               data = '\t'+hit.title
               if count <= self.hit_limit:
                   print >>data_file,data
                   gi, acc = hit.hit_id.split('|')[1:4:2]
                   self.refs[gi] = acc
                   print data
               if  self.DB:
                    self.load_sql(hsp_results)

    
    def InitDB(self):
        
        sql = self.DB_c.execute("""CREATE TABLE IF NOT EXISTS genotype_hits
            (align_length REAL,
             accession VARCHAR,
             hit_def VARCHAR,
             hit_id REAL,
             bits REAL,
             expect REAL,
             gaps REAL,
             identities INT,
             positives INT,
             seq_id VARCHAR,
             query_end INT,
             query_start INT,
             sbjct_end INT,
             sbjct_start INT,
             score REAL);""")
        
        self.DB_conn.commit()

            
    def load_sql(self, hsp_results):
        # {'accession': u'KF673377',
        #   'hit_def': u'Raphidocelis contorta strain SAG 11.81 18S ribosomal RNA gene, complete sequence',
        #   'hit_id': u'gi|664686838|gb|KF673377.1|',
        #   'hsps': [<Bio.Blast.Record.HSP object at 0x7f3175331d10>],
        #   'length': 2094,
        #   'title': u'gi|664686838|gb|KF673377.1| Raphidocelis contorta strain SAG 11.81 18S ribosomal RNA gene, complete sequence'}
        hsp_data = dict((k, str(v).translate(None, """"'(){}""")) for k,v in  hsp_results.items() )
                         
        sql ="""INSERT OR IGNORE INTO genotype_hits
            (align_length,
             accession,
             hit_def,
             hit_id,
             bits,
             expect,
             gaps,
             identities,
             positives,
             seq_id,
             query_end,
             query_start,
             sbjct_end,
             sbjct_start,
             score)
             values 
             ('{align_length}',
             '{accession}',
             '{hit_def}',
             '{hit_id}',
             '{bits}',
             '{expect}',
             '{gaps}',
             '{identities}',
             '{positives}',
             '{seq_id}',
             '{query_end}',
             '{query_start}',
             '{sbjct_end}',
             '{sbjct_start}',
             '{score}');""".format(**hsp_data)
        
        self.DB_c.execute(sql)



        
if __name__ == '__main__':
    blast = ParseBlast()
    if blast.DB:
        blast.InitDB()
    blast.get_data()
