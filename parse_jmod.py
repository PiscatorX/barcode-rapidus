#!/bin/python
import sys
import itertools as it

__author__ = "Andrew Ndhlovu"
__copyright__ = "Copyright 2018"
__license__ = "GPL"
__version__ = "3"
__maintainer__ = "Andrew Ndhlovu"
__email__ = "drewxdvst@outlook.com"


def init_param(jmod_out):

    fp  = open(jmod_out)

    data, pos = [ (line, fp.tell()) for line in fp if line.startswith('Arguments') ][0]

    seq_file_name = data.split()[3]

    return seq_file_name, fp, pos


def get_param(*args):

    start = False

    fp.seek(pos)

    for line in fp:
        if line.startswith('::Best Models::'):
            start = True
        if start:
            next(fp)
            params = ['Test']+ fp.next().split()
            next(fp)
            break

    return fp, params 


def get_values(fp, params, seq_file_name):
    print("#!/bin/bash")
    get_data = lambda x: list(zip(params, x.split()))
    update_vals = {'TrN+I+G':'TN93'}
    for data in map(get_data, fp):
        phyml_param = dict(data)
        if len(phyml_param) < 4:
            continue
        for k,v in phyml_param.items():
            for k2  in update_vals:
                if v == k2:
                   phyml_param[k] = update_vals[v]

        print('{Test}\tphyml -i {} --model {Model}  -f {f(a)},{f(c)},{f(g)},{f(t)}  -ts/tv {titv}  -b 1000 -v {pInv}  -a {gamma}  -s  BEST'.format(seq_file_name, **phyml_param))
    print("""
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
""")

seq_file_name, fp, pos = init_param(sys.argv[1])        
fp, params = get_param(fp, pos)
get_values(fp, params, seq_file_name

)
