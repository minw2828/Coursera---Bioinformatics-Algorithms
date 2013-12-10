#################################################################
#
# Author: Min Wang (san-heng-yi-shu@163.com)
#
# Date Created:
# 21 Oct 2013
#
# Coursera - Bioinformatics Algorithms
# - An Explosion of Hidden Messages
#  
# Clump Finding Problem: Find patterns forming clumps in a string.
#     Input: A string Genome, and integers k, L, and t.
#     Output: All distinct k-mers forming (L, t)-clumps in Genome.
#
# Example:
#
# Sample Input:
#      CGGACTCGACAGATGTGAAGAACGACAATGTGAAGACTCGACACGACAGAGTGAAGAGAAGAGGAAACATTGTAA
#      5 50 4
# 
# Sample Output:
#      CGACA GAAGA
#
####################################################################

import sys
import re
import itertools

def read_file(filename):
    f = open(filename, 'r')
    data = f.read()
    return data
    f.close()

def possible_kmers(base1,base2,k):
    bases = base1 + base2
    possible = set()
    for p in itertools.permutations(bases,k):
        if len(''.join(p))==9:
            possible.add(''.join(p))
    return possible

def frequency(seq,kmer,t):
    if seq.count(kmer) < t:
        return False
    return True

def filter_kmer(raw_kmers):
    fil_kmers = set()
    for kmer in raw_kmers:
        if frequency(seq,kmer,t) == True:
            fil_kmers.add(kmer)
    return fil_kmers

def kmer_position(seq,fil_kmer,t,L):
    positions = [m.start() for m in re.finditer(fil_kmer,seq)]
    if positions[t-1] - positions[0] < L:
        return True
    return False

def result(fil_kmers):
    result = set()
    for kmer in fil_kmers:
        if kmer_position(seq,kmer,t,L) == True:
            result.add(kmer)
    return result
        
if __name__ == '__main__':

    seq_filename = sys.argv[-3]
    seq = read_file(seq_filename).strip()
    base1_filename = sys.argv[-2]
    base1 = read_file(base1_filename).strip().split(' ')
    base2_filename = sys.argv[-1]
    base2 = read_file(base2_filename).strip().split(' ')
    k = 2
    L = 500
    t = 3
    # k,L,t = map(int,data[1].strip().split(' '))
    raw_kmers = possible_kmers(base1,base2, k)
    fil_kmers = filter_kmer(raw_kmers)
    final_result = result(fil_kmers)

    fw = open('./output.'+sys.argv[-4][:-3]+'.txt','w')
    fw.write(' '.join(final_result))
    fw.close()
