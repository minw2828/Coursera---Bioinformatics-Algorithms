#!/usr/bin/python

###############################################################################
#
# Author:
#
# Sanyk28 (san-heng-yi-shu@163.com)
#
# Date created:
#
# 10 Dec 2013
#
# Contig Generation Problem: 
# Generate the contigs from a collection of reads (with imperfect coverage).
#      Input: A collection of k-mers Patterns.
#      Output: All contigs in DeBruijn(Patterns).
#
# CODE CHALLENGE: Solve the Contig Generation Problem.
# 
# Sample Input:
#      ATG
#      ATG
#      TGT
#      TGG
#      CAT
#      GGA
#      GAT
#      AGA
#
# Sample Output:
#      AGA ATG ATG CAT GAT TGGA TGT
#
############################################################################## 

import sys
import timeit
import math
import re
import heapq
import random
import numpy as np
from itertools import combinations,product,izip,ifilter,chain
from collections import Counter,defaultdict
import debruijn_graph_from_kmers

def read_file(input_file):
    '''
    >>> pattern = read_file('test.contig_generation.txt')
    >>> pattern
    ['ATG', 'ATG', 'TGT', 'TGG', 'CAT', 'GGA', 'GAT', 'AGA']
    >>> pattern = read_file('test.contig_generation.extra.txt')
    '''
    f = open(input_file)
    data = f.readlines()
    f.close()
    return [item.strip() for item in data]

def build_de_bruijn_graph_from_kmers(pattern):
    '''
    >>> build_de_bruijn_graph_from_kmers(pattern)
    [('AG', ('GA',)), ('CA', ('AT',)), ('GG', ('GA',)), ('AT', ('TG',)), ('GA', ('AT',)), ('TG', ('GG', 'GT'))]
    '''
    return debruijn_graph_from_kmers.de_bruijn_graph(pattern).items()

def inout_degree(pattern):
    graph = build_de_bruijn_graph_from_kmers(pattern)
    outdegree = Counter([item[0] for item in graph])
    indegree = Counter([','.join(item[1]) for item in graph])
    keys = set(outdegree.keys()+indegree.keys())
    d = {}
    for key in keys:
        d[key] = (indegree[key],outdegree[key])
    return d

def prefix(kmer_pair):
    return [item[:-1] for item in kmer_pair]

def suffix(kmer_pair):
    return [item[1:] for item in kmer_pair]

def select(pattern):
    d = inout_degree(pattern)
    r = []
    for key,value in d.iteritems():
        value1,value2 = value
        if value1 != 1 or value2 > 1:
            r.append(key)
    return r

def path_graph(kmer_pairs):
    data = [(kmer_pair1,kmer_pair2) for kmer_pair1,kmer_pair2 in product(kmer_pairs,kmer_pairs) if suffix(kmer_pair1) == prefix(kmer_pair2)]
    edges = find_path(data)
    nodes = [prefix(item) for item in edges]+[suffix(edges[-1])]
    return (edges,nodes)

def find_ender(edges):
    item0 = [item[0] for item in edges]
    item1 = [item[1] for item in edges]
    c0 = Counter(item0)
    c1 = Counter(item1)
    result = [item for item in c1 if c0[item] < c1[item]]
    return random.choice(result)

def inverted_dictionary(edges):
    return dict([(item2,item1) for item1,item2 in edges])

def find_path(data):
    ivd = inverted_dictionary(data)
    ender = find_ender(data)
    path = [ender]
    while len(ivd) != 0:
        try:
            starter = ivd[ender]
            del ivd[ender]
            path.append(starter)
            ender = starter
        except:
            break
    return path[::-1]

def list_duplicates(data_list):
    data_list = [tuple(item) for item in data_list]
    tally = defaultdict(list)
    for i,item in enumerate(data_list):
        tally[item].append(i)
    return ((key,locs) for key,locs in tally.items() if len(locs) > 1)

def form_string(List):
    string = List[0] + ''.join([item[-1] for item in List[1:]])
    return string

def result(filename):
    pattern = read_file(filename)
    edges,nodes = path_graph(kmer_pairs)
    k = len(kmer_pairs[0][0])
    former = [item[0] for item in edges]
    latter = [item[1] for item in edges]
    results = form_string(former)[:k+d]+form_string(latter)
    return results

if __name__ == "__main__":

    start = timeit.default_timer()
    results = result(sys.argv[-1])
    fw = open('output.'+sys.argv[-1][:-4]+'.txt','w')
    fw.write(results)
    fw.close()
    stop = timeit.default_timer()
    print stop - start
