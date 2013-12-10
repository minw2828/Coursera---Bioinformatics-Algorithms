#################################################################
#
# Author: Min Wang (san-heng-yi-shu@163.com)
#
# Date Created:
# 13 Nov 2013
#
# Coursera - Bioinformatics Algorithms
# - Peculiar Statistics of the Forward and Reverse Half-Strands
#
# Minimum Skew Problem: 
#     Find a position in a genome minimizing the skew.
#     Input: A DNA string Genome.
#     Output: All integer(s) i minimizing Skew(Prefixi (Text)) 
#             among all values of i (from 0 to |Genome|).
#
# CODE CHALLENGE: Solve the Minimum Skew Problem.
#
# Sample Input:
#     TAAAGACTGCCGAGAGGCCAACACGAGTGCTAGAACGAGGGGCGTAAACGCGGGTCCGAT
#
# Sample Output:
#     11 24
#  
###################################################################

import sys

def read_file(filename):
    f = open(filename, 'r')
    genome = f.read()
    return genome
    f.close()

def skew(genome):
    return genome.count('G') - genome.count('C')

def prefix(genome,i):
    return genome[:i]

def skew_diagram(genome):
    skew_values = []
    for i in range(len(genome)+1):
        skew_values.append(skew(prefix(genome,i)))
    return skew_values

def minimum_skew(skew_values):
    result = []
    i = 0
    while i < len(skew_values):
        if skew_values[i] == min(skew_values):
            result.append(i)
        i += 1
    return result

if __name__ == '__main__':

    genome = read_file(sys.argv[-1]).strip().upper()
    skew_values = skew_diagram(genome)
    result = minimum_skew(skew_values)
    print ' '.join(map(str,result))

