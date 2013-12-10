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
# CODE CHALLENGE: Implement GIBBSSAMPLER.
#      Input: Integers k, t, and N, followed by a collection of strings Dna.
#      Output: The strings BestMotifs resulting from running GIBBSSAMPLER(Dna, k, t, N) 
#              with 20 random starts. Remember to use pseudocounts!
#
# Sample Input:
#      8 5 100
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
# Note: As with RANDOMIZEDMOTIFSEARCH, there is a very small chance that your 
#       algorithm may be implemented correctly but not return the correct answer. 
#       We suggest running your algorithm on another dataset if you do not get 
#       the correct answer the first time.
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
    >>> k,t,N,Dna = read_file('test.gibbssampler.txt')
    '''
    f = open(input_file)
    data = [item.strip() for item in f.readlines()]
    k,t,N = map(int,data[0].split(' '))
    f.close()
    return (k,t,N,data[1:])

def correct(dna,k):
    l = len(dna)-k+1
    return [dna[i:i+k] for i in range(l)]

def Kmers_array(k,Dna):
    return np.asarray([correct(dna,k) for dna in Dna])

def random_kmer_selection(k,t,l,kmers_array):
    return [kmers_array[i,random.randrange(l-k+1)] for i in range(t)]

def Profile(List,k,t):
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
    List = np.asarray([list(item) for item in List])
    i,score = 0,0
    while i < k:
        c = Counter(List[:,i])
        score += t - c.most_common(1)[0][1]
        i += 1
    return score

def compute(kmer,profile):
    l = len(kmer)
    order = ['A','C','G','T']
    return reduce(operator.mul,[profile.item(order.index(kmer[i]),i) for i in range(l)],1)

def pmpkp(kmer_array,k,profile):
    length = len(kmer_array)
    return heapq.nlargest(1,[(kmer_array[i],compute(kmer_array[i],profile)) for i in range(length)],key=lambda x:x[1])[0][0]

def Random(t):
    result = [random.random() for i in range(t)]
    s = sum(result)
    result = map(lambda x:x/s,result)
    l = len(result)
    xk = np.arange(l)
    custm = stats.rv_discrete(name='custm',values=(xk,result))
    R = custm.rvs(size=1)
    return R[0]

def prgkst(k,kmer_array,profile):
    '''
    Note: 'prgkst' is short for 'profile_randomly_generated_kmer_in_a_string_text'
    '''
    length = len(kmer_array)
    result = [(kmer_array[i],compute(kmer_array[i],profile)) for i in range(length)]
    probability_distribution = [item[1] for item in result]
    s = sum(probability_distribution)
    probability_distribution = map(lambda x:x/float(s),probability_distribution)
    l = len(probability_distribution)
    xk = np.arange(l)
    custm = stats.rv_discrete(name='custm',values=(xk,probability_distribution))
    R = custm.rvs(size=1)
    index = R[0]
    return result[index][0]

def gibbssampler(k,t,N,l,kmers_array):
    bestmotifs = random_kmer_selection(k,t,l,kmers_array)
    score_bestmotifs = Score(bestmotifs,k,t)
    motifs = random_kmer_selection(k,t,l,kmers_array)
    for j in range(N):
        i = Random(t)
        motifs.pop(i)
        profile = Profile(motifs,k,t)
        motifs_i = prgkst(k,kmers_array[i],profile)
        motifs.insert(i,motifs_i)
        score_motifs = Score(motifs,k,t)
        if score_motifs  < score_bestmotifs:
            bestmotifs = motifs
            score_bestmotifs = score_motifs
    return (bestmotifs,score_bestmotifs)

def result(filename):
    k,t,N,Dna = read_file(filename)
    if N > 1000:
        N = 1000
    l = len(Dna[0])
    kmers_array = Kmers_array(k,Dna)
    bm = []
    lowest_bestscore = float('inf')
    lowest_motifs = []
    i = 0
    while True:
        print 'i: '+str(i)
        print 'lowest_bestscore: '+str(lowest_bestscore)
        print 'lowest_motifs: '+str(lowest_motifs)
        bestmotifs,bestscore = gibbssampler(k,t,N,l,kmers_array)
        if bestscore < lowest_bestscore:
            lowest_bestscore = bestscore
            lowest_motifs = bestmotifs
            i = 0
        else:
            i += 1
        if i == 20:
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
