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
# CODE CHALLENGE: Solve the Overlap Graph Problem (restated below).
#      Input: A collection Patterns of k-mers.
#      Output: The overlap graph Overlap(Patterns), in the form of an adjacency list.
#
# Sample Input:
#      ATGCG
#      GCATG
#      CATGC
#      AGGCA
#      GGCAT
#
# Sample Output:
#      AGGCA -> GGCAT
#      CATGC -> ATGCG
#      GCATG -> CATGC
#      GGCAT -> GCATG
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
    >>> kmers = read_file('test.overlap_graph.txt')
    '''
    f = open(input_file)
    data = [item.strip() for item in f.readlines()]
    f.close()
    return data

def overlap_graph(data):
    n = len(data[0])-1
    return [(k1,k2) for k1,k2 in product(data,data) if k1 != k2 and k1.endswith(k2[:n])]

def result(filename):
    kmers = read_file(filename)
    results = overlap_graph(kmers)
    return results

if __name__ == "__main__":

    start = timeit.default_timer()
    results = result(sys.argv[-1])
    fw = open('output.'+sys.argv[-1][:-4]+'.txt','w')
    for r in results:
        fw.write('{0} -> {1}'.format(r[0],r[1])+'\n')
    fw.close()
    stop = timeit.default_timer()
    print stop - start
