#!/usr/bin/python

###############################################################################
#
# Author:
#
# Sanyk28 (san-heng-yi-shu@163.com)
#
# Date created:
#
# 28 Nov 2013
#
# Profile-most Probable k-mer Problem: 
# Find a Profile-most probable k-mer in a string.
#      Input: A string Text, an integer k, and a k * 4 matrix Profile.
#      Output: A Profile-most probable k-mer in Text.
#
# CODE CHALLENGE: Solve the Profile-most Probable k-mer Problem.
#
# Sample Input:
#      ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT
#      5
#      A C G T
#      0.2 0.4 0.3 0.1
#      0.2 0.3 0.3 0.2
#      0.3 0.1 0.5 0.1
#      0.2 0.5 0.2 0.1
#      0.3 0.1 0.4 0.2
#
# Sample Output:
#      CCGAG
#
############################################################################### 

import sys
import timeit
import regex
import heapq
import operator
import numpy as np
from itertools import combinations,product,izip,ifilter,chain
from collections import Counter,defaultdict

def read_file(input_file):
    f = open(input_file)
    data = [item.strip() for item in f.readlines()]
    f.close()
    return (data[0],int(data[1]),data[2].split(' '),np.asarray([map(float,item.split(' ')) for item in data[3:]]))

def correct(seq,k):
    return set(seq[i:i+k] for i in range(len(seq)-k+1))

def compute(kmer,order,profile):
    '''
    >>> kmer = 'GCCTA'
    >>> order = ['A', 'C', 'G', 'T']
    >>> profile
    array([[ 0.2,  0.4,  0.3,  0.1],
           [ 0.2,  0.3,  0.3,  0.2],
           [ 0.3,  0.1,  0.5,  0.1],
           [ 0.2,  0.5,  0.2,  0.1],
           [ 0.3,  0.1,  0.4,  0.2]])
    >>> compute(kmer,order,profile)
    0.00027
    '''
    c,i = [],0
    while i < len(kmer):
        c.append(profile.item(i,order.index(kmer[i])))
        i+=1
    return reduce(operator.mul,c,1)

def pmpkp(text,k,order,profile):
    corrects = correct(text,k)
    return [(kmer, compute(kmer,order,profile)) for kmer in corrects]

def result(filename):
    text,k,order,profile = read_file(filename)
    results = pmpkp(text,k,order,profile)
    return sorted(results,key=lambda x:x[1],reverse=True)[0][0]

if __name__ == "__main__":

    start = timeit.default_timer()
    results = result(sys.argv[-1])
    print results
    stop = timeit.default_timer()
    print stop - start
