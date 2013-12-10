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
# CODE CHALLENGE: Solve the String Composition Problem.
#      Input: An integer k and a string Text.
#      Output: Compositionk(Text), where the k-mers are written in lexicographic order.
#
# Sample Input:
#      5
#      CAATCCAAC
#
# Sample Output:
#      AATCC
#      ATCCA
#      CAATC
#      CCAAC
#      TCCAA
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
    >>> k,Text = read_file('test.string_composition.txt')
    '''
    f = open(input_file)
    data = [item.strip() for item in f.readlines()]
    f.close()
    return (int(data[0]),data[1])

def correct(dna,k):
    l = len(dna)-k+1
    return [dna[i:i+k] for i in range(l)]

def result(filename):
    k,Text = read_file(filename)
    results = sorted(correct(Text,k))
    return results

if __name__ == "__main__":

    start = timeit.default_timer()
    results = result(sys.argv[-1])
    print '\n'.join(results)
    fw = open('output.'+sys.argv[-1][:-4]+'.txt','w')
    fw.write('\n'.join(results))
    fw.close()
    stop = timeit.default_timer()
    print stop - start
