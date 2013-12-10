#!/usr/bin/python

###############################################################################
#
# Author:
#
# Sanyk28 (san-heng-yi-shu@163.com)
#
# Date created:
#
# 25 Nov 2013
#
# CODE CHALLENGE: Implement CYCLOPEPTIDESEQUENCING (pseudocode reproduced below).
#
# Sample Input:
#     0 113 128 186 241 299 314 427
#
# Sample Output:
#     186-128-113 186-113-128 128-186-113 128-113-186 113-186-128 113-128-186
#
############################################################################### 

import sys
import timeit
import operator
from itertools import permutations,product
from collections import Counter

table = {
    'G': 57,  'A': 71,  'S': 87,  'P': 97,  'V': 99,
    'T': 101, 'C': 103, 'I': 113, 'L': 113, 'N': 114, 
    'D': 115, 'K': 128, 'Q': 128, 'E': 129, 'M': 131, 
    'H': 137, 'F': 147, 'R': 156, 'Y': 163, 'W': 186}

def read_file(input_file):
    '''
    >>> Spectrum = read_file('test.cyclo_peptide_sequencing.extra.txt')
    >>> Spectrum
    [0, 71, 97, 99, 103, 113, 113, 114, 115, 131, 137, 196, 200, 202, 208, 214, 226, 227, 228, 240, 245, 299, 311, 311, 316, 327, 337, 339, 340, 341, 358, 408, 414, 
     424, 429, 436, 440, 442, 453, 455, 471, 507, 527, 537, 539, 542, 551, 554, 556, 566, 586, 622, 638, 640, 651, 653, 657, 664, 669, 679, 685, 735, 752, 753, 754, 
     756, 766, 777, 782, 782, 794, 848, 853, 865, 866, 867, 879, 885, 891, 893, 897, 956, 962, 978, 979, 980, 980, 990, 994, 996, 1022, 1093]
    '''
    f = open(input_file)
    raw_input = [int(item) for item in f.read().strip().split(' ')]
    f.close()
    return raw_input

# Stage 1
def select_single_peptides_spectrum(table,Spectrum):
    '''
    >>> Spectrum
    [0, 71, 97, 99, 103, 113, 113, 114, 115, 131, 137, 196, 200, 202, 208, 214, 226, 227, 228, 240, 245, 299, 311, 311, 316, 327, 337, 339, 340, 341, 358, 408, 414, 
     424, 429, 436, 440, 442, 453, 455, 471, 507, 527, 537, 539, 542, 551, 554, 556, 566, 586, 622, 638, 640, 651, 653, 657, 664, 669, 679, 685, 735, 752, 753, 754, 
     756, 766, 777, 782, 782, 794, 848, 853, 865, 866, 867, 879, 885, 891, 893, 897, 956, 962, 978, 979, 980, 980, 990, 994, 996, 1022, 1093]
    >>> sp = select_single_peptides_spectrum(table,Spectrum)
    >>> sp
    [71, 97, 99, 103, 113, 113, 114, 115, 131, 137]
    '''
    return [sp for sp in Spectrum if sp in table.values()]

# Stage 2
def initialize(sp,Spectrum):
    '''
    >>> Spectrum
    [0, 71, 97, 99, 103, 113, 113, 114, 115, 131, 137, 196, 200, 202, 208, 214, 226, 227, 228, 240, 245, 299, 311, 311, 316, 327, 337, 339, 340, 341, 358, 408, 414, 
     424, 429, 436, 440, 442, 453, 455, 471, 507, 527, 537, 539, 542, 551, 554, 556, 566, 586, 622, 638, 640, 651, 653, 657, 664, 669, 679, 685, 735, 752, 753, 754, 
     756, 766, 777, 782, 782, 794, 848, 853, 865, 866, 867, 879, 885, 891, 893, 897, 956, 962, 978, 979, 980, 980, 990, 994, 996, 1022, 1093]
    >>> sp = [71, 97, 99, 103, 113, 113, 114, 115, 131, 137]
    >>> initialize(sp,Spectrum)
    [(71, 131), (71, 137), (97, 99), (97, 103), (97, 131), (99, 97), (99, 103), (99, 115), (103, 97), (103, 99), (103, 137), (113, 113), (113, 113), (113, 114), 
     (113, 115), (113, 113), (113, 113), (113, 114), (113, 115), (114, 113), (114, 113), (114, 114), (114, 131), (115, 99), (115, 113), (115, 113), (131, 71), 
     (131, 97), (131, 114), (137, 71), (137, 103)]
    >>> len(initialize(sp,Spectrum))
    31
    '''
    return [t for t in product(sp,sp) if sum(t) in Spectrum]

# Stage 3
def stage3(subspectrum_list, sp, Spectrum):
    '''
    >>> subspectrum_list = initialize(sp,Spectrum)
    >>> subspectrum_list
    [(71, 131), (71, 137), (97, 99), (97, 103), (97, 131), (99, 97), (99, 103), (99, 115), (103, 97), (103, 99), (103, 137), (113, 113), (113, 113), (113, 114), 
     (113, 115), (113, 113), (113, 113), (113, 114), (113, 115), (114, 113), (114, 113), (114, 114), (114, 131), (115, 99), (115, 113), (115, 113), (131, 71), 
     (131, 97), (131, 114), (137, 71), (137, 103)]
    >>> sp
    [71, 97, 99, 103, 113, 113, 114, 115, 131, 137]
    >>> stage3(subspectrum_list, sp, Spectrum)
    [[71, 131, 97], [71, 131, 114], [71, 137, 103], [97, 99, 103], [97, 99, 115], [97, 103, 99], [97, 103, 137], [97, 131, 71], [99, 97, 103], [99, 97, 131], 
     [99, 103, 97], [99, 103, 137], [99, 115, 113], [103, 97, 99], [103, 99, 97], [103, 137, 71], [113, 113, 114], [113, 113, 115], [113, 114, 113], [113, 114, 131], 
     [113, 115, 99], [113, 115, 113], [114, 113, 113], [114, 131, 71], [115, 99, 97], [115, 113, 113], [131, 71, 137], [131, 97, 99], [131, 114, 113], [137, 71, 131], 
     [137, 103, 97], [137, 103, 99]]
    >>> len(stage3(subspectrum_list, sp, Spectrum))
    32
    '''
    poten = []
    for pt in product(subspectrum_list, sp):
        subspectrum = list(pt[0])+[pt[1]]
        skip = False
        c = Counter(subspectrum)
        for key,value in c.iteritems():
            if Spectrum.count(key) < value:
                skip = True
                break
        if skip == False and subspectrum_validity(subspectrum,Spectrum) != [] and subspectrum not in poten:
            poten.append(subspectrum)
    return poten

def subspectrum_validity(subspectrum,Spectrum):
    '''
    >>> Spectrum
    [0, 71, 97, 99, 103, 113, 113, 114, 115, 131, 137, 196, 200, 202, 208, 214, 226, 227, 228, 240, 245, 299, 311, 311, 316, 327, 337, 339, 340, 341, 358, 408, 414, 
     424, 429, 436, 440, 442, 453, 455, 471, 507, 527, 537, 539, 542, 551, 554, 556, 566, 586, 622, 638, 640, 651, 653, 657, 664, 669, 679, 685, 735, 752, 753, 754, 
     756, 766, 777, 782, 782, 794, 848, 853, 865, 866, 867, 879, 885, 891, 893, 897, 956, 962, 978, 979, 980, 980, 990, 994, 996, 1022, 1093]
    >>> subspectrum_validity([71, 131, 71])
    []
    >>> subspectrum_validity([103,137,71])
    [240, 208, 311]
    '''
    l = len(subspectrum)
    sv = [sum(subspectrum[start:start+2]) for start in range(0,l-1)] + [sum(subspectrum)]
    for item in sv:
        if item not in Spectrum:
            sv = []
            break
    return sv

# Stage 4
def iteration(sp,Spectrum):
    subspectrum_list = initialize(sp,Spectrum)
    i = 0
    while i < len(sp):
        subspectrum_list = stage3(subspectrum_list, sp, Spectrum)
        i = len(subspectrum_list[0])
    return subspectrum_list

# Stage 5
def generate_subpeptides(peptide,Spectrum):
    '''
    >>> generate_subpeptides([137, 71, 131, 114, 113, 113, 115, 99, 103, 97],Spectrum)
    [0, 71, 97, 99, 103, 113, 113, 114, 115, 131, 137, 200, 202, 202, 208, 214, 226, 227, 228, 234, 245, 299, 305, 316, 317, 327, 337, 339, 340, 341, 358, 408, 
     414, 429, 430, 436, 436, 440, 453, 455, 471, 507, 527, 539, 542, 543, 550, 551, 554, 566, 586, 622, 638, 640, 653, 657, 657, 663, 664, 679, 685, 735, 752, 
     753, 754, 756, 766, 776, 777, 788, 794, 848, 859, 865, 866, 867, 879, 885, 891, 891, 893, 956, 962, 978, 979, 980, 980, 990, 994, 996, 1022, 1093]
    '''
    subpeptides = [0,Spectrum[-1]]
    l = len(peptide)
    looped = peptide + peptide
    for start in range(0,l):
        for length in range(1,l):
            subpeptides.append(sum(looped[start:start+length]))
    return sorted(subpeptides)

def result(filename):
    Spectrum = read_file(filename)
    sp = select_single_peptides_spectrum(table,Spectrum)
    psc =  iteration(sp,Spectrum)
    results = []
    for item in psc:
        if generate_subpeptides(item,Spectrum) == Spectrum:
            results.append('-'.join(map(str,item)))
    return results

if __name__ == "__main__":

    start = timeit.default_timer()
    results = result(sys.argv[-1])
    fw = open('output.'+sys.argv[-1][:-4]+'.txt','w')
    fw.write(' '.join(results))
    fw.close()
    stop = timeit.default_timer()
    print stop - start
