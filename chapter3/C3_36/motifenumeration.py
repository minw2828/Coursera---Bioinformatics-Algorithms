#!/usr/bin/python

###############################################################################
#
# Author:
#
# Sanyk28 (san-heng-yi-shu@163.com)
#
# Date created:
#
# 27 Nov 2013
#
# CODE CHALLENGE: Implement MOTIFENUMERATION (reproduced below).
#      Input: Integers k and d, followed by a collection of strings Dna.
#      Output: All (k, d)-motifs in Dna.
#
# Sample Input:
#      3 1
#      ATTTGGC
#      TGCCTTA
#      CGGTATC
#      GAAAATT
#
# Sample Output:
#      ATA ATT GTT TTT
# 
# Note:
# Given a collection of strings Dna and an integer d, a k-mer is a (k,d)-motif 
# if it appears in every string from Dna with at most d mismatches.
#
############################################################################### 

import sys
import timeit
import regex
from itertools import combinations,product

def read_file(input_file):
    f = open(input_file)
    data = [item.strip() for item in f.readlines()]
    k,d = map(int,data[0].split(' '))
    f.close()
    return (k,d,data[1:])

def correct(seq,k):
    return set(seq[i:i+k] for i in range(len(seq)-k+1))

def correct_kmers(Dna,k):
    return frozenset().union(*[correct(seq,k) for seq in Dna])

def generate(s,d):
    N = len(s)
    letters = 'ACGT'
    pool = list(s)
    for indices in combinations(range(N),d):
        for replacements in product(letters,repeat=d):
            skip = False
            for i, a in zip(indices, replacements):
                if pool[i] == a:
                    skip = True
            if skip:
                continue
            key = dict(zip(indices,replacements))
            yield ''.join([pool[i] if i not in indices else key[i] for i in range(N)])

def possible_kmers(k,d,Dna):
    correct = set(correct_kmers(Dna,k))
    possibles = set()
    dd = 1
    while dd <= d:
        for s in correct:
            for item in generate(s,dd):
                possibles.add(item)
        dd += 1
    return possibles

def find_kmer(seq,kmer,d):
    return regex.findall(r'(?=(%s){s,e<=%d})'%(kmer,d),seq)

def kmer_composition(k,d,Dna):
    possibles = possible_kmers(k,d,Dna)
    kmers = []
    for kmer in possibles:
        skip = False
        for seq in Dna:
            if len(find_kmer(seq,kmer,d))<=0:
                skip = True
                break
        if skip == False:
            kmers.append(kmer)
    return kmers

def result(filename):
    k,d,Dna = read_file(filename)
    return kmer_composition(k,d,Dna)

if __name__ == "__main__":

    start = timeit.default_timer()
    results = result(sys.argv[-1])
    fw = open('output.'+sys.argv[-1][:-4]+'.txt','w')
    fw.write(' '.join(map(str,results)))
    fw.close()
    stop = timeit.default_timer()
    print stop - start
