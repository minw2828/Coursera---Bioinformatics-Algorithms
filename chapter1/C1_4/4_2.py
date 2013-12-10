##############################################################################
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
#################################################################################

import sys

def read_file(filename):
    f = open(filename, 'r')
    data = f.readlines()
    return data
    f.close()

def seq_window(seq, L):
    windows = set()
    for i in range(len(seq)-L+1):
        windows.add(seq[i:i+L])
    return windows

def possible_kmers(seq,k):
    possible = set()
    for i in range(len(seq)-k+1):
        possible.add(seq[i:i+k])
    return possible

def clump(seq,kmers,t):
    kmer = [item for item in kmers if seq.count(item) >= t]
    return kmer

def result(seq,k,L,t):
    kmer = []
    windows = seq_window(seq, L)
    for window in windows:
        kmers = possible_kmers(window,k)
        if len(clump(window,kmers,t)) > 0:
            kmer.append(' '.join(clump(window,kmers,t)))
    return kmer
        
if __name__ == '__main__':

    filename = sys.argv[-1]
    data  = read_file(filename)
    seq = data[0].strip()
    k = int(data[1].strip().split(' ')[0])
    L = int(data[1].strip().split(' ')[1])
    t = int(data[1].strip().split(' ')[2])
    
    raw_result =  ' '.join(result(seq,k,L,t))
    raw_result = raw_result.split(' ')
    final_result = set()
    for item in raw_result:
        final_result.add(item)
    
    fw = open('./output.'+filename,'w')
    fw.write(' '.join(final_result))
    fw.close()
