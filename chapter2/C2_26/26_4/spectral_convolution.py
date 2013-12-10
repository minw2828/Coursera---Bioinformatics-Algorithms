#!/usr/bin/python

###############################################################################
#
# Author:
#
# Sanyk28 (san-heng-yi-shu@163.com)
#
# Date created:
#
# 19 Nov 2013
#
# Spectral Convolution Problem: Compute the convolution of a spectrum.
#     Input: A collection of integers Spectrum.
#     Output: The list of elements in the convolution of Spectrum. If an element 
#             has multiplicity k, it should appearexactly k times; you may return 
#             the elements in any order.
#
# CODE CHALLENGE: Solve the Spectral Convolution Problem.
#
# Sample Input:
#     0 137 186 323
#
# Sample Output:
#     137 137 186 186 323 49
#
############################################################################### 

import sys
from itertools import product

def read_file(input_file):
    f = open(input_file)
    data = [int(item) for item in f.read().strip().split(' ')]
    f.close()
    return data

def result(filename):
    Spectrum = read_file(filename)
    return [pt1 - pt2 for pt1,pt2 in product(Spectrum,Spectrum) if pt1 != pt2 and pt1-pt2>=0]

if __name__ == "__main__":

    fw = open('output.'+sys.argv[-1][:-4]+'.txt','w')
    fw.write(' '.join(map(str,sorted(result(sys.argv[-1])))))
    fw.close()
