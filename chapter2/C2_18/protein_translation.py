#!/usr/bin/python

###############################################################################
#
# Author:
#
# Sanyk28 (san-heng-yi-shu@163.com)
#
# Date created:
#
# 18 Nov 2013
#
# Protein Translation Problem: Translate an RNA string into an amino acid string.
#     Input: An RNA string Pattern.
#     Output: The translation of Pattern into an amino acid string Peptide.
#
# CODE CHALLENGE: Solve the Protein Translation Problem.
#
# Sample Input:
#     AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA
#
# Sample Output:
#     MAMAPRTEINSTRING
################################################################################

import sys
from Bio.Seq import translate

def read_file(input_file):
    f = open(input_file)
    raw_input = f.read().strip()
    f.close()
    return raw_input

def translation(DNA):
    return translate(DNA, to_stop=True)

def result(filename):
    seq = read_file(filename)
    results = translation(seq)
    return results

if __name__ == '__main__':

    fw = open('output.txt','w')
    fw.write(result(sys.argv[-1]))
    fw.close()
