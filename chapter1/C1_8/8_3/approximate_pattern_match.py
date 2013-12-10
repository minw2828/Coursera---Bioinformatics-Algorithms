#######################################################################################
#
# Author: Min Wang (san-heng-yi-shu@163.com)
#
# Date Created:
# 13 Nov 2013
#
# Coursera - Bioinformatics Algorithms
# - Some Hidden Messages are More Elusive than Others
#
# Approximate Pattern Matching Problem: 
# Find all approximate occurrences of a pattern in a string.
#      Input: Two strings Pattern and Text along with an integer d.
#      Output: All positions where Pattern appears in Text with at most d mismatches.
#
# CODE CHALLENGE: Solve the Approximate Pattern Matching Problem
#
# Sample Input:
#     ATTCTGGA
#     CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAAT
#     3
#
# Sample Output:
#     6 7 26 27
#  
########################################################################################

import sys

def read_file(filename):
    f = open(filename, 'r')
    data = f.readlines()
    return data
    f.close()

def comparison(pattern,text):
    i,count = 0,0    
    while i < len(pattern):
        if pattern[i] != text[i]:
            count += 1
        i += 1
    return count

def approximate(pattern,text,d):
    i, result = 0, []
    while i < len(text)-len(pattern)+1:
        if comparison(pattern,text[i:i+len(pattern)]) <= d:
            result.append(i)
        i += 1
    return result

if __name__ == '__main__':

    pattern, text, d = [item.strip() for item in read_file(sys.argv[-1])]
    d = float(d)
    result = approximate(pattern,text,d)
    
    fw = open('output.'+sys.argv[-1][:-4]+'.txt','w')
    fw.write(' '.join(map(str,result)))
    fw.close()
