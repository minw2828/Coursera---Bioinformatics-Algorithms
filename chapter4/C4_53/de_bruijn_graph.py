#!/usr/bin/python

###############################################################################
#
# Author:
#
# Sanyk28 (san-heng-yi-shu@163.com)
#
# Date created:
#
# 3 Dec 2013
#
# CODE CHALLENGE: Solve the De Bruijn Graph from a String Problem.
#      Input: An integer k and a string Text.
#      Output: DeBruijnk(Text).
#
# Sample Input:
#      4
#      AAGATTCTCTAC
#
# Sample Output:
#      AAG -> AGA
#      AGA -> GAT
#      ATT -> TTC
#      CTA -> TAC
#      CTC -> TCT
#      GAT -> ATT
#      TCT -> CTA,CTC
#      TTC -> TCT
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
    >>> k,Text = read_file('test.de_bruijn_graph.txt')
    '''
    f = open(input_file)
    data = [item.strip() for item in f.readlines()]
    f.close()
    return (int(data[0]),data[1])

def correct(dna,k):
    l = len(dna)-k+1
    return [dna[i:i+k] for i in range(l)]

def overlap_graph(k,Text):
    data = correct(Text,k)
    n = len(data[0])-1
    return [(k1,k2) for k1,k2 in product(data,data) if k1 != k2 and k1.endswith(k2[:n])]

def de_bruijn_graph(k,Text):
    og = overlap_graph(k,Text)
    nodes = [Text[:k-1]]+[item[0][1:] for item in og]+[Text[-(k-1):]]
    raw_dbg = [(nodes[i],nodes[i+1]) for i in range(len(nodes)-1) if nodes[i] != nodes[i+1]]
    d = defaultdict(tuple)
    for tup in raw_dbg:
        d[tup[0]] += (tup[1],)
    for k,v in d.iteritems():
        if len(v) > 1:
            v = [(item,Text.index(item)) for item in v]
            v = sorted(v,key=lambda x:x[1],reverse=True)
            d[k] = [item[0] for item in v]
    return d
    
def result(filename):
    k,Text = read_file(filename)
    results = de_bruijn_graph(k,Text)
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
