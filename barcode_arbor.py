#! /usr/bin/env python
from Bio.Align.Applications import MuscleCommandline
from Bio.Alphabet import generic_dna
from Bio import AlignIO
import itertools as it
import subprocess
import argparse
import sys
import os




class BarcodeArbour(object):

    def __init__(self, args):
        
        self.threads = args.threads 
        self.msa_in = args.msa_in
        self.msa_out = args.msa_out if args.msa_out\
                       else ''.join([args.msa_in.replace('fasta',''),'aln'])

        self.jModelTest_jar = "java -jar /opt/jmodeltest2/dist/jModelTest.jar"
        #self.jModelTest_jar = "java -jar /home/andhlovu/bin/jmodeltest2/dist/jModelTest.jar"
        self.test_format = 'phylip'
        self.test_infile = ''.join([args.msa_in.replace('fasta',''), self.test_format ])
        self.test_out = ''.join([args.msa_in.replace('fasta',''), 'jModelTest' ])
        
    def run_muscle(self):
        
        # -in <inputfile>    Input file in FASTA format (default stdin)
        # -out <outputfile>  Output alignment in FASTA format (default stdout)
        # -diags             Find diagonals (faster for similar sequences)
        # -maxiters <n>      Maximum number of iterations (integer, default 16)
        # -maxhours <h>      Maximum time to iterate in hours (default no limit)
        # -html              Write output in HTML format (default FASTA)
        # -msf               Write output in GCG MSF format (default FASTA)
        # -clw               Write output in CLUSTALW format (default FASTA)
        # -clwstrict         As -clw, with 'CLUSTAL W (1.81)' header
        # -log[a] <logfile>  Log to file (append if -loga, overwrite if -log)
        # -quiet             Do not write progress messages to stderr
        # -version           Display version information and exit
        self.muscle_cline = MuscleCommandline(input=self.msa_in, out=self.msa_out)
        self.run_cmd(str(self.muscle_cline))
        assert os.path.exists(self.msa_out),'Did not find MSA file: {}'.format(self.msa_out)
        
        with open(self.msa_out) as input_handle, open(self.test_infile, "w") as output_handle:
            AlignIO.convert(input_handle, 'fasta', output_handle, self.test_format , generic_dna)    
        


    def run_jModelTest(self):

       # -d inputFile Sets the input data file. jModelTest makes use of the ALTER library for converting several alignment formats to PHYLIP
       # -s Sets the number of substitution schemes.
       # -f Include models with unequals base frecuencies
       # -i Include models with a proportion invariable sites.
       # -g numberOfRateCategories Include models with rate variation among sites and sets the number of categories. Usually 4 categories are enough.
       # -AIC Calculate the Akaike Information Criterion. 
       # -AICc Calculate the corrected Akaike Information Criterion.
       # -BIC Calculate the Bayesian Information Criterion. 
       # -DT Calculate the decision theory criterion. 
       # -p Calculate the parameter importances
       # -a Estimate model-averaged phylogeny for each active criterion.
       # -w Prints out the PAUP block.
       # -tr numberOfThreads Number of threads to execute (default is the number of logical processors in the machine).
       args = "-d {} -s 11 -f -i -g 4 -AIC -BIC -AICc -DT -p -a -w -tr {} -o {}".format(self.test_infile,\
                                                                                        self.threads,\
                                                                                        self.test_out)               
       jModelTest_cmd = ' '.join([self.jModelTest_jar, args])
       self.run_cmd(jModelTest_cmd)

    
       
    def parse_jModelTest(self):
    
       self.test_out = "AB1-TW4_final.jModelTest"
       fp = open(self.test_out)
       
       data, pos = [ (line, fp.tell()) for line in fp if line.startswith('Arguments') ][0]
       seq_file_name = data.split()[3]
       start = False

       fp.seek(pos)
       for line in fp:
           if line.startswith('::Best Models::'):
               start = True
           if start:
               fp.next()
               params = ['Test']+ fp.next().split()
               fp.next()
               break
       
       get_data = lambda x: zip(params, x.split())
       update_vals = {'TrN':'TN93'}
       phyl_cmd_fp = open(self.test_out+'.log','w')
       strip_model =  lambda model: model.split('+', 1)[0]
       del_val = {'titv': '-ts/tv',
                  'pInv': '-v',
                 'gamma': '-a'}
       
       self.model_data = {}
       cmd='{Test}\tphyml -i {} --model {Model} -f {f(a)},{f(c)},{f(g)},{f(t)} -ts/tv {titv} -b 1000 -v {pInv} -a {gamma} -s BEST'
       
       for data in map(get_data, fp):
           cmd_args = cmd
           del_values = []
           if len(data) < 4:
               continue
           phyml_param = dict(data)
           for k,v in phyml_param.items():
               if  v == 'N/A':
                   del_values = del_values + [k, del_val[k]]
               elif k == 'Model':
                   phyml_model = strip_model(phyml_param[k])
                   phyml_param[k] = update_vals.get(phyml_model, phyml_model)
           for val in del_values:
               cmd_args = cmd_args.replace(val, '')

           cmd_args = cmd_args.replace('{}','').format(seq_file_name, **phyml_param)
           self.model_data[phyml_param['Test']] = cmd_args.split("\t",1)[1]
           print >>phyl_cmd_fp, cmd_args
       log="""
#-i seq_file_name 
#-m (or --model) model
     # model : substitution model name.
     # - Nucleotide-based models : HKY85 (default) | JC69 | K80 | F81 | F84 | TN93 | GTR | custom (*)
     # (*) : for the custom option, a string of six digits identifies the model. For instance, 000000
     #  corresponds to F81 (or JC69 provided the distribution of nucleotide frequencies is uniform).
     # 012345 corresponds to GTR. This option can be used for encoding any model that is a nested within GTR.
#-f "fA,fC,fG,fT" : only valid for nucleotide-based models. fA, fC, fG and fT are floating numbers that 
     # correspond to the frequencies of A, C, G and T respectively (WARNING: do not use any blank space between
     # your values of nucleotide frequencies, only commas!)
#-t (or --ts/tv) ts/tv_ratio
     # ts/tv_ratio : transition/transversion ratio. DNA sequences only.
     # Can be a fixed positive value (ex:4.0) or e to get the maximum likelihood estimate.
#-v (or --pinv) prop_invar
     # prop_invar : proportion of invariable sites.
     # Can be a fixed value in the [0,1] range or e to get the maximum likelihood estimate.
# -a (or --alpha) gamma
     # gamma : distribution of the gamma distribution shape parameter.
     # Can be a fixed positive value or e to get the maximum likelihood estimate.
# -s (or --search) move
# 	  Tree topology search operation option.
# 	  Can be either NNI (default, fast) or SPR (a bit slower than NNI) or BEST (best of NNI and SPR search).
"""
       print >>phyl_cmd_fp, log
       print self.model_data['AICc']

       
       
    def run_cmd(self, cmd):
        print cmd
        process = subprocess.Popen(cmd,
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,
                                 shell=True)
        
        stdout, stderr= process.communicate()
        print stdout
        print stderr

        
    

if  __name__ ==  '__main__':
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("msa_in")
    parser.add_argument('-o','--msa-out')
    parser.add_argument('-t','--threads', type=int, default=2)
    args = parser.parse_args()
    arbour = BarcodeArbour(args)
    #arbour.run_muscle()
    #arbour.run_jModelTest()
    arbour.parse_jModelTest()
