#!/usr/bin/python

###############################################################################
#
# Author:
#
# Sanyk28 (san-heng-yi-shu@163.com)
#
# Date created:
#
# 29 Nov 2013
#
# CODE CHALLENGE: Implement GREEDYMOTIFSEARCH.
#
# Input: Integers k and t, followed by a collection of strings Dna.
#
# Output: A collection of strings BestMotifs resulting from applying 
#         GREEDYMOTIFSEARCH(Dna,k,t). If at any step you find more than one 
#         Profile-most probable k-mer in a given string, use the one occurring first.
#
# Sample Input:
#      3 5
#      GGCGTTCAGGCA
#      AAGAATCAGTCA
#      CAAGGAGTTCGC
#      CACGTCAATCAC
#      CAATAATATTCG
#
# Sample Output:
#      CAG
#      CAG
#      CAA
#      CAA
#      CAA
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
    '''
    >>> k,t,Dna = read_file('test.greedymotifsearch.txt')
    >>> Dna
    ['GGCGTTCAGGCA', 'AAGAATCAGTCA', 'CAAGGAGTTCGC', 'CACGTCAATCAC', 'CAATAATATTCG']
    >>> k,t,Dna = read_file('test.greedymotifsearch.extra.txt')
    '''
    f = open(input_file)
    data = [item.strip() for item in f.readlines()]
    k,t = map(int,data[0].split(' '))
    f.close()
    return (k,t,data[1:])

def correct(seq,k):
    return [seq[i:i+k] for i in range(len(seq)-k+1)]

def all_kmers(k,Dna):
    '''
    >>> Dna
    ['GGCGTTCAGGCA', 'AAGAATCAGTCA', 'CAAGGAGTTCGC', 'CACGTCAATCAC', 'CAATAATATTCG']
    >>> all_kmers(k,Dna)
    [['GGC', 'GCG', 'CGT', 'GTT', 'TTC', 'TCA', 'CAG', 'AGG', 'GGC', 'GCA'], 
     ['AAG', 'AGA', 'GAA', 'AAT', 'ATC', 'TCA', 'CAG', 'AGT', 'GTC', 'TCA'], 
     ['CAA', 'AAG', 'AGG', 'GGA', 'GAG', 'AGT', 'GTT', 'TTC', 'TCG', 'CGC'], 
     ['CAC', 'ACG', 'CGT', 'GTC', 'TCA', 'CAA', 'AAT', 'ATC', 'TCA', 'CAC'], 
     ['CAA', 'AAT', 'ATA', 'TAA', 'AAT', 'ATA', 'TAT', 'ATT', 'TTC', 'TCG']]
    '''
    return [correct(seq,k) for seq in Dna]

def profile_baseline_consensus(List,k,t):
    '''
    >>> List = ['GGC', 'AAG', 'CAA', 'CAC', 'CAA']
    >>> k = 3
    >>> t = 5
    >>> profile_baseline_consensus(List,k,t)
    ('CAA', 6, array([[ 0.2,  0.8,  0.4],
                      [ 0.6,  0. ,  0.4],
                      [ 0. ,  0. ,  0. ],
                      [ 0.2,  0.2,  0.2]]))
    >>> List = ['GATGGACCGGGC', 'TAAAAAGGTATA', 'AACCACGAGTAC', 'TGTCATGTGCGG', 'AACCTAAACCCT', 'GTGCCCGATATG', 'TAGTCTTCGAGG', 'AGGAGACGTGTT', 'TGTGGGGATCGT', 'CGCAGTGCACTA', 'TACTCGTAACTT', 'GCTCTAGTACGC', 'GAGACGGTCGTA', 'GATCGGTGGCAG', 'GTAGGTATCACC', 'GTGGCTATCGCT', 'TGAGCAGACCCG', 'AGTGATCTGAGC', 'CAAAATGGGAGT', 'GTTGGTATCACC', 'CCTCGGAAAACG', 'GGCGGCTCCATC', 'TACTAGTATAAG', 'GTGGTTATCACC', 'CATCACGCAATG']
    >>> k = 12
    >>> t = 25
    >>> profile_baseline_consensus(List,k,t)
    ('CGTGATCACTGA', 190, array([[ 0.16,  0.28,  0.16,  0.12,  0.48,  0.24,  0.08,  0.32,  0.24, 0.16,  0.16,  0.28],
                                 [ 0.4 ,  0.08,  0.28,  0.28,  0.04,  0.32,  0.4 ,  0.12,  0.28, 0.32,  0.28,  0.28],
                                 [ 0.12,  0.16,  0.32,  0.24,  0.24,  0.36,  0.28,  0.24,  0.28, 0.36,  0.2 ,  0.2 ],
                                 [ 0.32,  0.48,  0.24,  0.36,  0.24,  0.08,  0.24,  0.32,  0.2 , 0.16,  0.36,  0.24]]))
    '''
    List = np.asarray([list(item) for item in List])
    i,baseline_score,consensus = 0,0,''
    pro = np.zeros(shape=(4,k))
    while i < k:
        d = {'A':0, 'C':0, 'G':0, 'T':0}
        c = Counter(List[:,i])
        consensus += c.most_common(1)[0][0]
        baseline_score += t-c.most_common(1)[0][1]
        for key,value in c.iteritems():
            d[key] = value/float(t)
        pro[:,i] = [d['A'],d['C'],d['G'],d['T']]
        i += 1
    return (consensus,baseline_score,pro)

def pmpkp(kmer,order,profile):
    '''
    >>> kmer = 'GGC'
    >>> order = ['A', 'C', 'G', 'T']
    >>> profile
    array([[ 0.2,  0.8,  0.4],
           [ 0.6,  0. ,  0.4],
           [ 0. ,  0. ,  0. ],
           [ 0.2,  0.2,  0.2]])
    >>> pmpkp(kmer,order,profile)
    0.0
    >>> kmer = 'CAT'
    >>> pmpkp(kmer,order,profile)
    0.096
    >>> kmer = 'GATGGACCGGGC'
    >>> profile
    array([[ 0.16,  0.28,  0.16,  0.12,  0.48,  0.24,  0.08,  0.32,  0.24, 0.16,  0.16,  0.28],
           [ 0.4 ,  0.08,  0.28,  0.28,  0.04,  0.32,  0.4 ,  0.12,  0.28, 0.32,  0.28,  0.28],
           [ 0.12,  0.16,  0.32,  0.24,  0.24,  0.36,  0.28,  0.24,  0.28, 0.36,  0.2 ,  0.2 ],
           [ 0.32,  0.48,  0.24,  0.36,  0.24,  0.08,  0.24,  0.32,  0.2 , 0.16,  0.36,  0.24]])
    >>> pmpkp(kmer,order,profile)
    3.0204666209894405e-08
    '''
    c,i = [],0
    while i < len(kmer):
        c.append(profile.item(order.index(kmer[i]),i))
        i+=1
    return reduce(operator.mul,c,1)

def form_profile(iks,kmer_set,k,t):
    '''
    >>> iks = ['GGC']
    >>> kmer_set = ['AAG', 'AGA', 'GAA', 'AAT', 'ATC', 'TCA', 'CAG', 'AGT', 'GTC', 'TCA']
    >>> form_profile(iks,kmer_set,k,t)
    ['GGC', 'AAG']
    '''
    order = ['A','C','G','T']
    consensus,baseline,profile = profile_baseline_consensus(iks,k,t)
    intermediate = [(kmer,pmpkp(kmer,order,profile)) for kmer in kmer_set]
    iks.append(heapq.nlargest(1,intermediate,key=lambda x:x[1])[0][0])
    return iks

def loopover(iks,ksl,k,t):
    '''
    >>> iks = ['GGC']
    >>> ksl = all_kmers(k,Dna)[1:]
    >>> iks = ['GTACATCTCTCT']
    >>> loopover(iks,ksl,k,t)
    ['GGC', 'AAG', 'TTC', 'GTC', 'TTC']
    >>> ksl = all_kmers(k,Dna)[1:]
    >>> loopover(iks,ksl,k,t)
    ['GTACATCTCTCT', 'GATAGTTGCCGT', 'CGCCAAGATTAC', 'CGGATACAGTGG', 'TTTAAGTACCAC', 'CGCAGGCGTGAT', 'GGTAGATATCCC', 
     'CTCAGACTGGTG', 'TTGATTGTGGCG', 'CTGAGGCGTTTG', 'CTTCAACATTAT', 'TTGATATGTTAT', 'GGTAGTGACTAG', 'CTCAAGTATCCC', 
     'GGTATGTGGCGC', 'GGTAGTCGGGAT', 'GTTAGTGGGTAT', 'CGGCTACTGTTC', 'TTCCTTCATGCC', 'CTCAGGTGTTAC', 'TGTAGGTATCAC', 
     'GGTCGTTATCCC', 'GTCATTCTGGGC', 'CGGCTTGATCTC', 'GGTAATCGGGAC']
    '''
    i = 0
    while i < len(ksl):
        iks = form_profile(iks,ksl[i],k,t)
        i +=1
    return iks

def best_motifs_collections(k,t,Dna):
    '''
    >>> Dna
    ['GGCGTTCAGGCA', 'AAGAATCAGTCA', 'CAAGGAGTTCGC', 'CACGTCAATCAC', 'CAATAATATTCG']
    >>> aksv = all_kmers(k,Dna)
    >>> aksv
    [['GGC', 'GCG', 'CGT', 'GTT', 'TTC', 'TCA', 'CAG', 'AGG', 'GGC', 'GCA'], 
     ['AAG', 'AGA', 'GAA', 'AAT', 'ATC', 'TCA', 'CAG', 'AGT', 'GTC', 'TCA'], 
     ['CAA', 'AAG', 'AGG', 'GGA', 'GAG', 'AGT', 'GTT', 'TTC', 'TCG', 'CGC'], 
     ['CAC', 'ACG', 'CGT', 'GTC', 'TCA', 'CAA', 'AAT', 'ATC', 'TCA', 'CAC'], 
     ['CAA', 'AAT', 'ATA', 'TAA', 'AAT', 'ATA', 'TAT', 'ATT', 'TTC', 'TCG']]
    >>> best_motifs_collections(k,t,Dna)
    [['GGC', 'AAG', 'TTC', 'GTC', 'TTC'], ['GCG', 'AAG', 'CAA', 'AAT', 'AAT'], ['CGT', 'AAG', 'AAG', 'AAT', 'AAT'], 
     ['GTT', 'AAG', 'AAG', 'AAT', 'AAT'], ['TTC', 'AAG', 'AGT', 'AAT', 'AAT'], ['TCA', 'AAG', 'CAA', 'CAA', 'CAA'], 
     ['CAG', 'AAG', 'CAA', 'CAA', 'CAA'], ['AGG', 'AAG', 'CAA', 'AAT', 'AAT'], ['GGC', 'AAG', 'TTC', 'GTC', 'TTC'], 
     ['GCA', 'TCA', 'CAA', 'TCA', 'CAA']]
    '''
    aksv = all_kmers(k,Dna)
    return [loopover([iks],aksv[1:],k,t) for iks in aksv[0]]

def greedy_motif_search(Dna,k,t):
    '''
    >>> Dna = ['GGCGTTCAGGCA', 'AAGAATCAGTCA', 'CAAGGAGTTCGC', 'CACGTCAATCAC', 'CAATAATATTCG']
    >>> k = 3
    >>> t = 5
    >>> aksv = all_kmers(k,Dna)
    >>> aksv
    [['GGC', 'GCG', 'CGT', 'GTT', 'TTC', 'TCA', 'CAG', 'AGG', 'GGC', 'GCA'], 
     ['AAG', 'AGA', 'GAA', 'AAT', 'ATC', 'TCA', 'CAG', 'AGT', 'GTC', 'TCA'], 
     ['CAA', 'AAG', 'AGG', 'GGA', 'GAG', 'AGT', 'GTT', 'TTC', 'TCG', 'CGC'], 
     ['CAC', 'ACG', 'CGT', 'GTC', 'TCA', 'CAA', 'AAT', 'ATC', 'TCA', 'CAC'], 
     ['CAA', 'AAT', 'ATA', 'TAA', 'AAT', 'ATA', 'TAT', 'ATT', 'TTC', 'TCG']]
    >>> best_motifs = [dna[:k] for dna in Dna]
    >>> best_motifs
    ['GGC', 'AAG', 'CAA', 'CAC', 'CAA']
    >>> greedy_motif_search(Dna,k,t)
    ['CAG', 'AAG', 'CAA', 'CAA', 'CAA']
    '''
    aksv = all_kmers(k,Dna)
    best_motifs = [dna[:k] for dna in Dna]
    baseline_score = profile_baseline_consensus(best_motifs,k,t)[1]
    bmcs = best_motifs_collections(k,t,Dna)
    for bmc in bmcs:
        pbc = profile_baseline_consensus(bmc,k,t)
        if pbc[1] < baseline_score:
            baseline_score = pbc[1]
            best_motifs = bmc
    return best_motifs

def result(filename):
    k,t,Dna = read_file(filename)
    return greedy_motif_search(Dna,k,t)

if __name__ == "__main__":

    start = timeit.default_timer()
    results = result(sys.argv[-1])
    fw = open('output.'+sys.argv[-1][:-4]+'.txt','w')
    fw.write('\n'.join(results))
    fw.close()
    stop = timeit.default_timer()
    print stop - start
