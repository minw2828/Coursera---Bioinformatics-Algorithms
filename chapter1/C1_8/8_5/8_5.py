###############################################################################
#
# Author: Min Wang (san-heng-yi-shu@163.com)
#
# Date Created:
# 21 Nov 2013
#
# Coursera - Bioinformatics Algorithms
# - Some Hidden Messages are More Elusive than Others
#
# Frequent Words with Mismatches and Reverse Complements Problem: 
# Find the most frequent k-mers (with mismatches and reverse complements) in a DNA string.
#       Input: A DNA string Text as well as integers k and d.
#       Output: All k-mers Pattern maximizing the sum Countd(Text, Pattern) + Countd(Text, Pattern)
#       over all possible k-mers.
#
# CODE CHALLENGE: 
# Solve the Frequent Words with Mismatches and Reverse Complements Problem.
#
# Sample Input:
#      ACGTTGCATGTCGCATGATGCATGAGAGCT
#      4 1
#
# Sample Output:
#      ATGT ACAT
#
###############################################################################

import sys
import regex
import timeit
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from itertools import combinations, product

def read_file(filename):
    f = open(filename, 'r')
    data = f.readlines()
    f.close()
    return data

def reverse_complement(seq):
    my_dna = Seq(seq, generic_dna)
    rc = my_dna.reverse_complement()
    return str(rc)

def correct_kmers(seq,k):
    correct = set(seq[i:i+k] for i in range(len(seq)-k+1))
    return correct

def generate_mismatches(s,d):
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
            for item in generate_mismatches(s,dd):
                possibles.add(item)
        dd += 1
    return possibles

def all_kmers_pairup(seq,k,d):
    possibles = possible_kmers(seq,k,d)
    all_pair_kmers = set()
    for kmer in possibles:
        rc_kmer = reverse_complement(kmer)
        if (kmer,rc_kmer) not in all_pair_kmers and (rc_kmer,kmer) not in all_pair_kmers:
            all_pair_kmers.add((kmer,reverse_complement(kmer)))
    return all_pair_kmers

def find_kmer(seq,kmer,d):
    match = regex.findall(r'(?=(%s){s,e<=%d})'%(kmer,d),seq)
    return len(match)

def kmer_composition(seq,k,d):
    all_pair_kmers = all_kmers_pairup(seq,k,d)
    new_dict = {}
    for pair_kmers in all_pair_kmers:
        new_dict[pair_kmers] = sum([find_kmer(seq,i,d) for i in pair_kmers])
    return new_dict

def result(filename):
    seq,num = [item.strip() for item in read_file(filename)]
    k,d = [int(item.strip()) for item in num.split(' ')]
    new_dict = kmer_composition(seq,k,d)
    maximum = max([value for value in new_dict.itervalues()])
    results = [' '.join(key) for key in new_dict.iterkeys() if new_dict[key] == maximum]
    return results

if __name__ == '__main__':

    start = timeit.default_timer()
    results = result(sys.argv[-1])
    stop = timeit.default_timer()
    print stop - start
    fw = open('output.'+sys.argv[-1][:-4]+'.txt','w')
    fw.write(' '.join(map(str,results)))
    fw.close()
