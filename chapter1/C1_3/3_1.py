#################################################################
#
# Author: Min Wang (san-heng-yi-shu@163.com)
#
# Date Created:
# 21 Oct 2013
#
# Coursera - Bioinformatics Algorithms
# - Some Hidden Messages are More Surprising than Others
# 
# Reverse Complement Problem: Reverse complement a nucleotide pattern.
#      Input: A DNA string Pattern.
#      Output: Pattern, the reverse complement of Pattern.
#
# Example:
#
# Sample Input:
#      AAAACCCGGT
# 
# Sample Output:
#      ACCGGGTTTT 
####################################################################

import sys

complement = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}

def read_file(filename):
    f = open(filename, 'r')
    data = f.readline().strip()
    return data
    f.close()

def reverse_complement(seq):
    seqc = [complement[ch] for ch in seq]
    seqc = ''.join(seqc)
    return seqc[::-1]

if __name__ == '__main__':

    argv = sys.argv[-1]
    seq = read_file(argv)

    fw = open('./output.3_1.txt','w')
    fw.write(reverse_complement(seq))
    fw.close()
