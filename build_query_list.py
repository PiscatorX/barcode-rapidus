#! /usr/bin/env python

import  argparse

__author__ = "Andrew Ndhlovu"
__copyright__ = "Copyright 2018"
__license__ = "GPL"
__version__ = "3"
__maintainer__ = "Andrew Ndhlovu"
__email__ = "drewxdvst@outlook.com"


def build_Q(q_list, query_builder):
    for q  in set(q_list):
        print(query_builder.replace("{}",q))


if __name__ ==  '__main__':        
    parser = argparse.ArgumentParser(description="build a GenBank query for each item on provided list e.g.  \"(Prorocentrum[Primary Organism] AND 28S[Text Word])\" Prorocentrum is the query and the query builder would be \"({}[Primary Organism] AND 28S[Text Word])\". ")
    parser.add_argument('query_list', action='store', type=str, help="list of search string to be inserted in query")
    parser.add_argument('-b', '--build-query', dest='query_builder', action='store',required=True, type=str, help="query builder use {} for item replacement in query term ")
    args  =  parser.parse_args()
    q_list =  open(args.query_list).read().splitlines()
    q = args.query_builder
    build_Q(q_list, q)
