#!/usr/bin/python

###############################################################################
#
# Author:
#
# Sanyk28 (san-heng-yi-shu@163.com)
#
# Date created:
#
# 20 Nov 2013
#
# CODE CHALLENGE: Implement CONVOLUTIONCYCLOPEPTIDESEQUENCING.
#
# Input: An integer M, an integer N, and a collection of (possibly repeated) 
#        integers Spectrum.
#
# Output: A cyclic peptide LeaderPeptide with amino acids taken only from the 
#         top M elements (and ties) of the convolution of Spectrum that fall 
#         between 57 and 200, and where the size of Leaderboard is restricted 
#         to the top N (and ties).
#
# Sample Input:
#     20
#     60
#     57 57 71 99 129 137 170 186 194 208 228 265 285 299 307 323 356 364 394 422 493
#
# Sample Output:
#     99-71-137-57-72-57
#
############################################################################### 

import sys

def read_file(input_file):
    f = open(input_file)
    raw_input = f.read()
    f.close()
    return raw_input

def result(filename):
    data = [int(item) for item in read_file(filename).strip().split(' ')]
    results = []
    for item1 in data:
        for item2 in data:
            if item1 != item2 and item1-item2 >= 0:
                results.append(item1-item2)
    return results
  
if __name__ == "__main__":

    fw = open('output.'+sys.argv[-1][:-4]+'.txt','w')
    fw.write(' '.join(map(str,sorted(result(sys.argv[-1])))))
    fw.close()
