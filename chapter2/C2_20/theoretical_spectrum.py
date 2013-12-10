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
# Generating Theoretical Spectrum Problem: 
# Generate the theoretical spectrum of a cyclic peptide.
#     Input: An amino acid string Peptide.
#     Output: Cyclospectrum(Peptide).
#
# CODE CHALLENGE: Solve the Generating Theoretical Spectrum Problem.
#
# Sample Input:
#     LEQN
#
# Sample Output:
#     0 113 114 128 129 227 242 242 257 355 356 370 371 484
#
############################################################################### 

import sys

def read_file(input_file):
    f = open(input_file)
    raw_input = f.readlines()
    f.close()
    return raw_input

def lookup_table(filename):
    table = {}
    for line in read_file(filename):
        table[line.strip().split(' ')[0]] = int(line.strip().split(' ')[1])
    return table

def generate_subpeptides(peptide):
    subpeptides = ['',peptide]
    l = len(peptide)
    looped = peptide + peptide
    for start in range(0,l):
        for length in range(1,l):
            subpeptides.append(looped[start:start+length])
    return subpeptides

def result(table_filename, peptide_filename):
    protein_mass_table = lookup_table(table_filename)
    peptide = read_file(peptide_filename)[0].strip()
    subpeptides = generate_subpeptides(peptide)
    results = []
    for item in subpeptides:
        try:
            mass = sum([protein_mass_table[c] for c in item])
            results.append(mass)
        except:
            results.append(0)
    return sorted(results)

if __name__ == "__main__":

    fw = open('output.'+sys.argv[-1][:-4]+'.txt','w')
    fw.write(' '.join(map(str,result('integer_mass_table.txt', sys.argv[-1]))))
    fw.close()
