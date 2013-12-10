#!/usr/bin/python

###############################################################################
#
# Author:
#
# Sanyk28 (san-heng-yi-shu@163.com)
#
# Date created:
#
# 2 Dec 2013
#
# CODE CHALLENGE: Implement RANDOMIZEDMOTIFSEARCH.
#      Input: Integers k and t, followed by a collection of strings Dna.
#      Output: A collection BestMotifs resulting from running 
#              RANDOMIZEDMOTIFSEARCH(Dna, k, t) 1000 times.
#              Remember to use pseudocounts!
#
# Sample Input:
#      8 5
#      CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA
#      GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG
#      TAGTACCGAGACCGAAAGAAGTATACAGGCGT
#      TAGATCAAGTTTCAGGTGCACGTCGGTGAACC
#      AATCCACCAGCTCCACGTGCAATGTTGGCCTA
#
# Sample Output:
#      TCTCGGGG
#      CCAAGGTG
#      TACAGGCG
#      TTCAGGTG
#      TCCACGTG
#
############################################################################### 

import sys
import timeit
import heapq
import operator
import random
import numpy as np
from itertools import combinations,product,izip,ifilter,chain
from collections import Counter,defaultdict

def read_file(input_file):
    '''
    >>> k,t,Dna = read_file('test.randomizedmotifsearch.txt')
    >>> k,t,Dna = read_file('test.randomizedmotifsearch.extra.txt')
    '''
    f = open(input_file)
    data = [item.strip() for item in f.readlines()]
    k,t = map(int,data[0].split(' '))
    f.close()
    return (k,t,data[1:])

def random_kmer_selection(k,t,Dna):
    '''
    >>> random_kmer_selection(k,t,Dna)
    ['CGGGGGTG', 'CGAGGTAT', 'AGACCGAA', 'CAGGTGCA', 'ACCAGCTC']
    >>> random_kmer_selection(k,t,Dna)
    ['CCCCTCTC', 'TGCCAAGG', 'TACAGGCG', 'CGTCGGTG', 'ATCCACCA']
    '''
    l = len(Dna[0])
    kmers = []
    for dna in Dna:
        n = random.randrange(l-k)
        kmers.append(dna[n:n+k])
    return kmers

def Profile(List,k,t):
    '''
    >>> List = ['CGGGGGTG', 'CGAGGTAT', 'AGACCGAA', 'CAGGTGCA', 'ACCAGCTC']
    >>> Profile(List,k,t)
    array([[ 0.5       ,  0.33333333,  0.5       ,  0.33333333,  0.16666667, 0.16666667,  0.5       ,  0.5       ],
           [ 0.66666667,  0.33333333,  0.33333333,  0.33333333,  0.33333333, 0.33333333,  0.33333333,  0.33333333],
           [ 0.16666667,  0.66666667,  0.5       ,  0.66666667,  0.66666667, 0.66666667,  0.16666667,  0.33333333],
           [ 0.16666667,  0.16666667,  0.16666667,  0.16666667,  0.33333333, 0.33333333,  0.5       ,  0.33333333]])
    '''
    List = np.asarray([list(item) for item in List])
    i = 0
    pro = np.ones(shape=(4,k))
    while i < k:
        c = Counter(List[:,i])
        pro[0,i] = pro[0,i]+c['A']
        pro[1,i] = pro[1,i]+c['C']
        pro[2,i] = pro[2,i]+c['G']
        pro[3,i] = pro[3,i]+c['T']
        i += 1
    return pro/float(t+1)

def Score(List,k,t):
    '''
    >>> List = ['CGGGGGTG', 'CGAGGTAT', 'AGACCGAA', 'CAGGTGCA', 'ACCAGCTC']
    >>> Score(List,k,t)
    19
    '''
    List = np.asarray([list(item) for item in List])
    i,score = 0,0
    while i < k:
        c = Counter(List[:,i])
        score += t - c.most_common(1)[0][1]
        i += 1
    return score

def compute(kmer,order,profile):
    c,i = [],0
    while i < len(kmer):
        c.append(profile.item(order.index(kmer[i]),i))
        i+=1
    return reduce(operator.mul,c,1)

def correct(seq,k):
    '''
    >>> seq = 'CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA'
    >>> correct(seq,k)
    set(['GTGTTCAG', 'CCTCTCGG', 'CAGTAAAC', 'CCCTCTCG', 'TCTCGGGG', 'TGTTCAGT', 'GGGTGTTC', 'GCCCCTCT', 'GTAAACGG', 
         'GGGGGTGT', 'GGGGTGTT', 'CGGGGGTG', 'GGTGTTCA', 'CGCCCCTC', 'AACGGCCA', 'GTTCAGTA', 'AGTAAACG', 'TCGGGGGT', 
         'TAAACGGC', 'TCAGTAAA', 'TTCAGTAA', 'CTCTCGGG', 'CCCCTCTC', 'CTCGGGGG', 'AAACGGCC'])
    '''
    return set([seq[i:i+k] for i in range(len(seq)-k+1)])

def pmpkp(dna,k,order,profile):
    '''
    >>> dna = 'AATCCACCAGCTCCACGTGCAATGTTGGCCTA'
    >>> pmpkp(dna,k,order,profile)
    'CACGTGCA'
    '''
    corrects = correct(dna,k)
    return heapq.nlargest(1,[(kmer,compute(kmer,order,profile)) for kmer in corrects],key=lambda x:x[1])[0][0]

def Motifs(profile,k,t,Dna):
    '''
    >>> motifs = ['CGGGGGTG', 'CGAGGTAT', 'AGACCGAA', 'CAGGTGCA', 'ACCAGCTC']
    >>> Profile(motifs,k,t)
    array([[ 0.5       ,  0.33333333,  0.5       ,  0.33333333,  0.16666667, 0.16666667,  0.5       ,  0.5       ],
           [ 0.66666667,  0.33333333,  0.33333333,  0.33333333,  0.33333333, 0.33333333,  0.33333333,  0.33333333],
           [ 0.16666667,  0.66666667,  0.5       ,  0.66666667,  0.66666667, 0.66666667,  0.16666667,  0.33333333],
           [ 0.16666667,  0.16666667,  0.16666667,  0.16666667,  0.33333333, 0.33333333,  0.5       ,  0.33333333]])
    >>> Motifs(profile,k,t,Dna)
    ['CGGGGGTG', 'CGAGGTAT', 'AGACCGAA', 'CAGGTGCA', 'CACGTGCA']
    '''
    order = ['A', 'C', 'G', 'T']
    return [pmpkp(dna,k,order,profile) for dna in Dna]

def random_motif_search(k,t,Dna):
    '''
    >>> list(random_motif_search(k,t,Dna))
    [['CAGTAAAC', 'AAGGTGCC', 'AAGAAGTA', 'CAGGTGCA', 'CACGTGCA']]
    '''
    bestmotifs = random_kmer_selection(k,t,Dna)
    initial_motifs = random_kmer_selection(k,t,Dna)
    while True:
        profile = Profile(initial_motifs,k,t)
        motifs = Motifs(profile,k,t,Dna)
        if Score(motifs,k,t) < Score(bestmotifs,k,t):
            bestmotifs = motifs
            initial_motifs = motifs
            # initial_motifs = random_kmer_selection(k,t,Dna)
        else:
            yield bestmotifs
            break

def result(filename):
    k,t,Dna = read_file(filename)
    bm = []
    lowest_bestscore = float('inf')
    lowest_motifs = []
    i = 0
    while True:
        bestmotifs = list(random_motif_search(k,t,Dna))[0]
        bestscore = Score(bestmotifs,k,t)
        if bestscore < lowest_bestscore:
            lowest_bestscore = bestscore
            lowest_motifs = bestmotifs
            i = 0
        else:
            i += 1
        if i > 900:
            break
    return lowest_motifs 
    
if __name__ == "__main__":

    start = timeit.default_timer()
    results = result(sys.argv[-1])
    print '\n'.join(results)
    fw = open('output.'+sys.argv[-1][:-4]+'.txt','w')
    fw.write('\n'.join(results))
    fw.close()
    stop = timeit.default_timer()
    print stop - start
