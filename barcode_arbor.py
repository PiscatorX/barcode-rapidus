#! /usr/bin/env python
from Bio.Align.Applications import MuscleCommandline
from Bio.Alphabet import generic_dna
from Bio import AlignIO
import itertools as it
import subprocess
import argparse
import pprint
import sys
import os




class BarcodeArbour(object):

    def __init__(self, args):
        
        self.threads = args.threads 
        self.msa_in = args.msa_in
        self.msa_out = args.msa_out if args.msa_out\
                       else ''.join([args.msa_in.replace('fasta',''),'aln'])

        self.jModelTest_jar = "java -jar /home/andhlovu/bin/jmodeltest2/dist/jModelTest.jar"
        self.modeltest_fmt = 'phylip'
        self.mrbayes_fmt =  'nexus'
        self.modeltest_infile = ''.join([args.msa_in.replace('fasta',''), self.modeltest_fmt ])
        self.modeltest_out = ''.join([args.msa_in.replace('fasta',''), 'jModelTest' ])
        self.mrbayes_infile = ''.join([args.msa_in.replace('fasta',''), self.mrbayes_fmt ])
        self.paup_blocks = {}    
        self.Test = {}
        self.model_args = {}
        self.phyml_template ='phyml -i {msa_fname} --model {Model} -f {f(a)},{f(c)},{f(g)},{f(t)} -ts/tv {titv} -b 1000 -v {pInv} -a {gamma} -s BEST --no-memory-check'
        self.mrbayes_params = {'nruns': 1,
                                'ngen': 100,
                          'samplefreq': 10,
                                'file': self.mrbayes_infile}

        
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
        self.conv(self.msa_out, self.modeltest_infile, 'fasta', self.modeltest_fmt)



        
    def conv(self, infile, outfile, informat, outformat):
        
        with open(infile) as input_handle, open(outfile, "w") as output_handle:
            AlignIO.convert(input_handle, informat, output_handle, outformat, generic_dna)    
        


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
       args = "-d {} -s 11 -f -i -g 4 -AICc -p -a -w -tr {} -o {}".format(self.modeltest_infile,\
                                                                                        self.threads,\
                                                                                        self.modeltest_out)
       
       jModelTest_cmd = ' '.join([self.jModelTest_jar, args])
       self.run_cmd(jModelTest_cmd)
    
    

    def parse_jModelTest(self):

       modelTest_fp = open(self.modeltest_out)       
       for line in modelTest_fp:
            self.block = []
            if line.startswith('Arguments'):
                self.seq_file_name = line.split()[3]
            if line.startswith("[!"):
                while not line.startswith('Lset'):
                     line = modelTest_fp.next().strip()
                     self.block.append(line)
                self.PauP2MrBayes()
            if line.startswith("::Best"):
                self.fp_pointer = modelTest_fp 
                self.GetModel()
       self.PhymlArgs()
       modelTest_fp.close()


       
    def GetModel(self):
        
        
        header = ['Test']+[ self.fp_pointer.next().strip() for i in range(2) ][1].split()
        line = self.fp_pointer.next()
        while any(line):
            line = self.fp_pointer.next().strip().split()
            if any(line):
                model_params = dict(zip(header, line))
                model_params.update({'msa_fname': self.seq_file_name})
                self.Test[model_params['Test']] = model_params
                 
                      
    
    def PhymlArgs(self):
        
        args  = {'titv': '-ts/tv',
                   'pInv': '-v',
                  'gamma': '-a'}
        model_name = {'TrN':'TN93'}
        get_model =  lambda model: model.split('+', 1)[0]
        phyl_cmd_fp = open(self.modeltest_out+'.log','w')
        for test_name, params_data in self.Test.items():
            cmd = self.phyml_template 
            rm_list = []
            for param, value in params_data.items():    
                if value == 'N/A':
                    params_data.update({param:''})
                    cmd = cmd.replace(args.get(param,param),'')
                if  param == 'Model':
                    model = get_model(value)
                    params_data['Model'] = model_name.get(model,model)
            cmd = cmd.format(**params_data)
            self.model_args[test_name] = cmd
            cmd_args = "{}\t{}".format(test_name, cmd)
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

        
            
    def PauP2MrBayes(self):

        model = self.block[0].split()[-1]
        block = [ item for word in self.block[4].split('=')\
                   for item in word.rsplit(" ",1) ][1:]
        data_block =  dict([(param,value) for param, value in\
                                    it.izip_longest(*[iter(block)]*2)])

        self.paup_blocks[model] = data_block
        self.mrbayes_params.update(data_block)
        mrbayes_cmd = """begin mrbayes;
        set autoclose=yes nowarn=yes;
        execute {file};
        lset nst={nst} rates={rates};
        mcmc nruns={nruns} ngen={ngen} samplefreq={samplefreq} file={file}1;
        mcmc file={file}2;
        mcmc file={file}3;
        end;""".format(**self.mrbayes_params)
        print mrbayes_cmd
        exit(1)


        
        # #NST (number of substitution types)
        # mrbayes_cmd = """begin mrbayes;
        # set autoclose=yes nowarn=yes;
        # execute {file};
        # lset nst={nst} rates={rates};
        # mcmc nruns={nruns} ngen={ngen} samplefreq={samplefre} file={file}1;
        # mcmc file={file}2;
        # mcmc file={file}3;
        # end;""".format(self.mrbayes_params)
        # 
        
        

        

        
        
        
    def phyml(self):

        phyml_cmd = self.model_args['AICc']
        self.run_cmd(phyml_cmd)

        
    def mrbayes(self):
        self.conv(self.msa_out, self.mrbayes_infile, 'fasta', self.mrbayes_fmt)



        

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
    arbour.run_muscle()
    arbour.run_jModelTest()
    arbour.parse_jModelTest()
    arbour.phyml()
    arbour.mrbayes()
