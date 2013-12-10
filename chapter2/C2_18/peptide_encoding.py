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
# Peptide Encoding Problem: 
# Find substrings of a genome encoding a given amino acid sequence.
#     Input: A DNA string Text and an amino acid string Peptide.
#     Output: All substrings of Text encoding Peptide (if any such substrings exist).
#
# CODE CHALLENGE: Solve the Peptide Encoding Problem.
#
# Sample Input:
#     ATGGCCATGGCCCCCAGAACTGAGATCAATAGTACCCGTATTAACGGGTGA
#     MA
#
# Sample Output:
#     ATGGCC
#     GGCCAT
#     ATGGCC
#
# Note: The solution may contain repeated strings if the same string occurs 
#       more than once as a substring of Text and encodes Peptide.
#
############################################################################### 

import sys
import re
from itertools import product
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

dna_codon_table = {
    'TTT': 'F',      'CTT': 'L',      'ATT': 'I',      'GTT': 'V',
    'TTC': 'F',      'CTC': 'L',      'ATC': 'I',      'GTC': 'V',
    'TTA': 'L',      'CTA': 'L',      'ATA': 'I',      'GTA': 'V',
    'TTG': 'L',      'CTG': 'L',      'ATG': 'M',      'GTG': 'V',
    'TCT': 'S',      'CCT': 'P',      'ACT': 'T',      'GCT': 'A',
    'TCC': 'S',      'CCC': 'P',      'ACC': 'T',      'GCC': 'A',
    'TCA': 'S',      'CCA': 'P',      'ACA': 'T',      'GCA': 'A',
    'TCG': 'S',      'CCG': 'P',      'ACG': 'T',      'GCG': 'A',
    'TAT': 'Y',      'CAT': 'H',      'AAT': 'N',      'GAT': 'D',
    'TAC': 'Y',      'CAC': 'H',      'AAC': 'N',      'GAC': 'D',
    'TAA': 'Stop',   'CAA': 'Q',      'AAA': 'K',      'GAA': 'E',
    'TAG': 'Stop',   'CAG': 'Q',      'AAG': 'K',      'GAG': 'E',
    'TGT': 'C',      'CGT': 'R',      'AGT': 'S',      'GGT': 'G',
    'TGC': 'C',      'CGC': 'R',      'AGC': 'S',      'GGC': 'G',
    'TGA': 'Stop',   'CGA': 'R',      'AGA': 'R',      'GGA': 'G',
    'TGG': 'W',      'CGG': 'R',      'AGG': 'R',      'GGG': 'G' }

def read_file(input_file):
    f = open(input_file)
    raw_input = f.readlines()
    f.close()
    return raw_input

def reverse_complement(dna):
    my_dna = Seq(dna,generic_dna)
    my_dna_rc = my_dna.reverse_complement()
    return str(my_dna_rc)

def generate_kmers(dna,peptide):
    codons, kmers = [],[]
    for d,v in dna_codon_table.iteritems():
        if v == peptide[0]:
            codons.append(d)

    for codon in codons:
        for index in [m.start() for m in re.finditer(r'(?=(%s))'%codon,dna)]:
            kmers.append(dna[index:index+len(peptide)*3])
    return kmers

def translation(DNA):
    my_dna = Seq(DNA, generic_dna)
    protein = my_dna.translate()
    return str(protein)

def confirm_tranlation(dna,peptide):
    kmers = generate_kmers(dna,peptide)
    confirmed_kmers = []
    for k in kmers:
        if translation(k) == peptide:
            confirmed_kmers.append(k)
    return confirmed_kmers

def result(filename):
    dna,peptide = [item.strip() for item in read_file(filename)]   
    results = [item for item in confirm_tranlation(dna,peptide)]
    for item in confirm_tranlation(reverse_complement(dna),peptide):
        results.append(reverse_complement(item))
    return results

if __name__ == "__main__":

    print '\n'.join(result(sys.argv[-1]))
