#!/usr/bin/python

###############################################################################
#
# Author:
#
# Sanyk28 (san-heng-yi-shu@163.com)
#
# Date created:
#
# 4 Dec 2013
#
# DeBruijn Graph from k-mers Problem: Construct the de Bruijn graph from a set of k-mers.
#      Input: A collection of k-mers Patterns.
#      Output: The adjacency list of the de Bruijn graph DeBruijn(Patterns).
#
# CODE CHALLENGE: Solve the de Bruijn Graph from k-mers Problem.
#
# Sample Input:
#      GAGG
#      GGGG
#      GGGA
#      CAGG
#      AGGG
#      GGAG
#
# Sample Output:
#      AGG -> GGG
#      CAG -> AGG
#      GAG -> AGG
#      GGA -> GAG
#      GGG -> GGA,GGG
#
############################################################################### 

import sys
import timeit
import heapq
import operator
import random
import numpy as np
from scipy import stats
from itertools import combinations,product,izip,ifilter,chain
from collections import Counter,defaultdict

def read_file(input_file):
    '''
    >>> kmers = read_file('test.debruijn_graph_from_kmers.txt')
    >>> kmers = read_file('test.debruijn_graph_from_kmers.extra.txt')
    '''
    f = open(input_file)
    data = [item.strip() for item in f.readlines()]
    f.close()
    return data

def overlap_graph(kmers):
    return [(k1,k2) for k1,k2 in product(kmers,kmers) if k1 != k2 and k1.endswith(k2[:-1])]

def de_bruijn_graph(kmers):
    og = overlap_graph(kmers)
    nodes = [(item[0][:-1],item[1][:-1]) for item in og]
    potentials = [(item[0][1:],item[1][1:]) for item in og]
    for item in potentials:
        if item not in nodes:
            nodes.append(item)
    d = defaultdict(tuple)
    for tup in nodes:
        if tup[1] not in d[tup[0]]:
            d[tup[0]] += (tup[1],)
    return d
    
def result(filename):
    kmers = read_file(filename)
    results = de_bruijn_graph(kmers)
    return results

if __name__ == "__main__":

    start = timeit.default_timer()
    results = result(sys.argv[-1])
    fw = open('output.'+sys.argv[-1][:-4]+'.txt','w')
    for k,v in results.iteritems():
        fw.write('{0} -> {1}'.format(k,','.join(v))+'\n')
    fw.close()
    stop = timeit.default_timer()
    print stop - start
