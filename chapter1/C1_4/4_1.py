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
# Pattern Matching Problem: 
#      Find all occurrences of a pattern in a string.
#      Input: Two strings, Pattern and Text.
#      Output: All starting positions where Pattern appears as a 
#              substring of Text.
#
# Example:
#
# Sample Input:
#      ATAT
#      GATATATGCATATACTT
# 
# Sample Output:
#      1 3 9
#
#################################################################

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

    filename, pattern = sys.argv[-2], sys.argv[-1]
    genome = read_file(filename)
    result = ' '.join(map(str,occurrences(pattern, genome)))
    
    fw = open('./output.4_1.'+pattern+'.txt', 'w')
    fw.write(result)
    fw.close()
