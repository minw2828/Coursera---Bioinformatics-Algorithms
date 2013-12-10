#################################################################
#
# Author: Min Wang (san-heng-yi-shu@163.com)
#
# Date Created:
# 21 Nov 2013
#
# Coursera - Bioinformatics Algorithms
# - Some Hidden Messages are More Elusive than Others
#
# Frequent Words with Mismatches Problem: 
#     Find the most frequent k-mers with mismatches in a string.
#     Input: A string Text as well as integers k and d. (You may 
#            assume k <= 12 and d <= 3.)
#     Output: All most frequent k-mers with up to d mismatches in Text.
#
# CODE CHALLENGE: Solve the Frequent Words with Mismatches Problem.
#
# Sample Input:
#     ACGTTGCATGTCGCATGATGCATGAGAGCT 4 1
# Sample Output:
#     GATG ATGC ATGT
#  
###################################################################

import sys
import regex
import timeit
from itertools import combinations, product

def read_file(filename):
    f = open(filename, 'r')
    data = f.read()
    f.close()
    return data

def correct_kmers(seq,k):
    correct = set(seq[i:i+k] for i in range(len(seq)-k+1))
    return correct

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

def possible_kmers(seq,k,d):
    possibles = set()
    correct = correct_kmers(seq,k)
    dd = 1
    while dd <= d:
        for s in correct:
            for item in generate(s,dd):
                possibles.add(item)
        dd += 1
    return possibles

def find_kmer(seq,kmer,d):
    match = regex.findall(r'(?=(%s){s,e<=%d})'%(kmer,d),seq)
    return match

def kmer_composition(seq,k,d):
    possibles = possible_kmers(seq,k,d)
    kmers = {}
    for kmer in possibles:
        kmers[kmer] = len(find_kmer(seq,kmer,d))
    return kmers

def result(filename):
    seq,k,d = [item.strip() for item in read_file(filename).split(' ')]
    k, d = int(k), int(d)
    kmers = kmer_composition(seq,k,d)
    maximum = max([value for value in kmers.itervalues()])
    results = [key for key in kmers.iterkeys() if kmers[key] == maximum]
    return results

if __name__ == '__main__':

    start = timeit.default_timer()    
    results = result(sys.argv[-1])
    stop = timeit.default_timer()
    print stop - start
    fw = open('output.'+sys.argv[-1][:-4]+'.txt','w')
    fw.write(' '.join(map(str,results)))
    fw.close()
