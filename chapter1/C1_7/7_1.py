###################################################################
#
# Author: Min Wang (san-heng-yi-shu@163.com)
#
# Date Created:
# 13 Nov 2013
#
# Coursera - Bioinformatics Algorithms
# - Peculiar Statistics of the Forward and Reverse Half-Strands
#
# EXERCISE BREAK: Give all values of Skew(Prefixi (GAGCCACCGCGATA)) 
#                 for i ranging from 0 to 14.
#
# Sample Input:
#     CATGGGCATCGGCCATACGCC
#
# Sample Output:
#     0 -1 -1 -1 0 1 2 1 1 1 0 1 2 1 0 0 0 0 -1 0 -1 -2
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
    result = []
    for i in range(15):
    # for i in range(len(genome)+1):
        result.append(skew(prefix(genome,i)))
    return result

if __name__ == '__main__':

    # genome = read_file(sys.argv[-1]).strip().upper()
    genome = 'GAGCCACCGCGATA'
    result = skew_diagram(genome)
    print ' '.join(map(str,result))

