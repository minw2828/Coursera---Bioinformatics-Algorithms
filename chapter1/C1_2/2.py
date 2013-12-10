#################################################################
#
# Author: Min Wang (san-heng-yi-shu@163.com)
#
# Date Created:
# 21 Oct 2013
#
# Coursera - Bioinformatics Algorithms
# - Hidden Messages in the Replication Origin
# 
# Frequent Words Problem: Find the most frequent k-mers in a string.
# Input: A string Text and an integer k.
# Output: All most frequent k-mers in Text.
#
# Example:
#
# Sample Input:
#      ACGTTGCATGTCGCATGATGCATGAGAGCT
#      4
#
# Sample Output:
#     CATG GCAT
# 
####################################################################

import sys

def read_file(filename):
    f = open(filename, 'r')
    data = f.readlines()
    return data
    f.close()

def possible_kmers(seq,n):
    possible = set()
    for i in range(len(seq)-n+1):
        possible.add(seq[i:i+n])
    return possible

def most_frequent_kmers(seq,kmers):
    freqs = [seq.count(kmer) for kmer in kmers]
    maximum = max(freqs)
    kmer = [kmer for kmer in kmers if seq.count(kmer) == maximum]
    return kmer

if __name__ == '__main__':

    argv = sys.argv[-1]
    seq, n = read_file(argv)[0].strip(), int(read_file(argv)[1].strip())
    kmers = possible_kmers(seq,n)
    result = most_frequent_kmers(seq,kmers)
    print ' '.join(result)
