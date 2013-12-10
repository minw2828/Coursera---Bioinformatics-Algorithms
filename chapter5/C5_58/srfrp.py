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
# CODE CHALLENGE: Solve the String Reconstruction from Read-Pairs Problem.
#      Input: An integer d followed by a collection of paired k-mers PairedReads.
#      Output: A string Text with (k, d)-mer composition equal to PairedReads.
#
# Sample Input:
#      2
#      GAGA|TTGA
#      TCGT|GATG
#      CGTG|ATGT
#      TGGT|TGAG
#      GTGA|TGTT
#      GTGG|GTGA
#      TGAG|GTTG
#      GGTC|GAGA
#      GTCG|AGAT
#
# Sample Output:
#      GTGGTCGTGAGATGTTGA
#
############################################################################## 

import sys
import timeit
import math
import random
from itertools import combinations,product,izip,ifilter,chain
from collections import Counter,defaultdict

def seperate_pairwise(pairs):
    return [tuple(item.split('|')) for item in pairs]

def read_file(input_file):
    '''
    >>> d,kmer_pairs = read_file('test.srfrp.txt')
    >>> d,kmer_pairs = read_file('test.srfrp.extra.txt')
    '''
    f = open(input_file)
    data = f.readlines()
    f.close()
    data = [item.strip() for item in data]
    return (int(data[0]),seperate_pairwise(data[1:]))

def prefix(kmer_pair):
    return [item[:-1] for item in kmer_pair]

def suffix(kmer_pair):
    return [item[1:] for item in kmer_pair]

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

def form_string(List):
    string = List[0] + ''.join([item[-1] for item in List[1:]])
    return string

def result(filename):
    d,kmer_pairs = read_file(filename)
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
