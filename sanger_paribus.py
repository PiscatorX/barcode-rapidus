#! /usr/bin/env python
from collections import Counter
import itertools as it
import argparse
import glob
import shutil
import os
import pprint

__author__ = "Andrew Ndhlovu"
__copyright__ = "Copyright 2018"
__license__ = "GPL"
__version__ = "3"
__maintainer__ = "Andrew Ndhlovu"
__email__ = "drewxdvst@outlook.com"
Entrez.email = "drewxdvst@outlook.com"



class  SangerParibus(object):

    def __init__(self):
        parser = argparse.ArgumentParser(description="""rename foward and reverse reads to match and create renamed copies for downstream merging""")
        parser.add_argument('-f', '--fwd', action='store', required=True, help='forward read directory')
        parser.add_argument('-r', '--rev', action='store', help='reverse read directory')
        parser.add_argument('-o', '--ouput_dir', action='store',help ='output directory')
        parser.add_argument('-t', '--tag', action='store',required=True, help ='uniques read tag for the filename')

        self.args = parser.parse_args()
        if not self.args.ouput_dir:
            self.args.ouput_dir = self.args.tag  

            
    def get_reads(self):
    
        self.fwd =  map(self.find_reads, [os.path.join(self.args.fwd, '*.ab1')])[0]
        self.rev =  map(self.find_reads, [os.path.join(self.args.rev, '*.ab1')])[0]
        self.uniq_ids = Counter(map(lambda x: x[0], self.fwd+self.rev))  

        
    def find_reads(self, file_patt):

        cut = lambda fname : (os.path.basename(fname.split('_')[0]), fname)
        return [ cut(read) for read in glob.iglob(file_patt) ] 

    
    def paribus(self):

        fetch = lambda read_id, read_list : [ read_ref[1] \
                for  read_ref in read_list if  read_id in read_ref ]

        for read_id in self.uniq_ids:
            fwd_fnames = fetch(read_id,  self.fwd)
            rev_fnames = fetch(read_id,  self.rev)
            i = 1
            print
            for fwd_read,  rev_read in it.izip_longest(fwd_fnames, rev_fnames):
                fwd_fname =  '{}_{}_{}_F.ab1'.format(read_id,self.args.tag,i)
                rev_fname =  '{}_{}_{}_R.ab1'.format(read_id,self.args.tag,i)
                map(self.mv,  [ (fwd_read, fwd_fname), (rev_read, rev_fname)])
                i += 1

                
    def mv(self,*args):
        
        read, fname  =  args[0]
        if read:
            dest_fname =os.path.join(self.args.ouput_dir, fname)
            dest_dir = self.args.ouput_dir
            try:
               os.makedirs(dest_dir)
            except OSError as e: 
               if os.path.isdir(dest_dir):
                  pass
               else:
                   raise Exception,e
            shutil.copy(read, dest_fname)
            

if  __name__ ==  '__main__':
    sanger_reads = SangerParibus()
    sanger_reads.get_reads()
    sanger_reads.paribus()
    
