#!/usr/bin/python

###############################################################################
#
# Author:
#
# Sanyk28 (san-heng-yi-shu@163.com)
#
# Date created:
#
# 5 Dec 2013
#
# CODE CHALLENGE: Solve the String Reconstruction Problem.
#      Input: The adjacency list of a directed graph that has an Eulerian path.
#      Output: An Eulerian path in this graph.
#
# Sample Input:
#      CTT -> TTA
#      ACC -> CCA
#      TAC -> ACC
#      GGC -> GCT
#      GCT -> CTT
#      TTA -> TAC
# 
# Sample Output:
#      GGCTTACCA
#
############################################################################### 

import sys
import timeit
import re
import heapq
import random
import numpy as np
from itertools import combinations,product,izip,ifilter,chain
from collections import Counter,defaultdict

def read_file(input_file):
    '''
    >>> data = read_file('test.string_reconstruction.txt')
    >>> data = read_file('test.string_reconstruction.extra.txt')
    '''
    f = open(input_file)
    data = [item.strip().split(' -> ') for item in f.readlines()]
    f.close()
    return data

def find_ender(data):
    item0 = [item[0] for item in data]
    item1 = [item[1] for item in data]
    c0 = Counter(item0)
    c1 = Counter(item1)
    result = [(item,c0[item],c1[item]) for item in c1 if c0[item] != c1[item]]
    return [item[0] for item in result if item[1] < item[2]][0]

def find_path(data):
    initial_ender = find_ender(data)
    ender = find_ender(data)
    path = [ender]
    while len(data) != 0:
        try:
            starter = [item1 for [item1,item2] in data if item2 == ender][0]
            data.remove([starter,ender])
            path.append(starter)
            ender = starter
        except:
            break
    return path[::-1]

def form_string(path):
    string = path[0]
    for item in path[1:]:
        string += item[-1]
    return string

def result(filename):
    data = read_file(filename)
    path = find_path(data)
    results = form_string(path)
    return results

if __name__ == "__main__":

    start = timeit.default_timer()
    results = result(sys.argv[-1])
    fw = open('output.'+sys.argv[-1][:-4]+'.txt','w')
    fw.write(results)
    fw.close()
    stop = timeit.default_timer()
    print stop - start
