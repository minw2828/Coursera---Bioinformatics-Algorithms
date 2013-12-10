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
# Pattern Matching Problem: Find all occurrences of a pattern in a string.
#      Input: Two strings, Pattern and Text.
#      Output: All starting positions where Pattern appears as a substring of Text.
#
# Example:
#
# Sample Input:
#      ATAT
#      GATATATGCATATACTT
# 
# Sample Output:
#      1 3 9
###################################################################

import sys
import re

def read_file(filename):
    f = open(filename, 'r')
    genome = f.readline().strip()
    return genome
    f.close()

def occurrences(pattern, genome):
    matches = re.finditer(r'(?=(%s))' % re.escape(pattern), genome)
    return [m.start(1) for m in matches]

if __name__ == '__main__':

    argv = sys.argv[-1]
    genome = read_file(argv)
    result = ' '.join(map(str,occurrences('CTTGATCAT', genome)))
    
    fw = open('./output.3_3.txt', 'w')
    fw.write(result)
    fw.close()
