#!/usr/bin/python

###############################################################################
#
# Author:
#
# Sanyk28 (san-heng-yi-shu@163.com)
#
# Date created:
#
# 27 Nov 2013
#
# CODE CHALLENGE: Implement CONVOLUTIONCYCLOPEPTIDESEQUENCING.
#
# Input: An integer M, an integer N, and a collection of (possibly repeated) 
#        integers Spectrum.
#
# Output: A cyclic peptide LeaderPeptide with amino acids taken only from the 
#         top M elements (and ties) of the convolution of Spectrum that fall 
#         between 57 and 200, and where the size of Leaderboard is restricted 
#         to the top N (and ties).
#
# Sample Input:
#      20
#      60
#      57 57 71 99 129 137 170 186 194 208 228 265 285 299 307 323 356 364 394 422 493
#
# Sample Output:
#      99-71-137-57-72-57
#
############################################################################### 

import sys
import timeit
import heapq
from itertools import product,izip,ifilter
from collections import Counter

def read_file(input_file):
    f = open(input_file)
    M,N,Spectrum = [item.strip() for item in f.readlines()]
    f.close()
    return (int(M),int(N),map(int,Spectrum.split(' ')))

def select_spectrum(Spectrum,M):
    poten = [pt1-pt2 for pt1,pt2 in product(Spectrum,Spectrum) if pt1-pt2 >= 57 and pt1-pt2 <= 200]
    lst = Counter(poten).most_common()
    tie = heapq.nlargest(M,lst,key=lambda x:x[1])[-1][1]
    res = list(ifilter(lambda x: x[1]>=tie, lst))
    return list(izip(*res))[0]

def score(theorectical_spectrum,experimental_spectrum):
    return len(list((Counter(theorectical_spectrum) & Counter(experimental_spectrum)).elements()))

def generate_subspectrums(peptide):
    l = len(peptide)
    looped = peptide + peptide
    return [0,sum(peptide)]+[sum(looped[start:start+length]) for start,length in product(range(0,l),range(1,l))]

def cut(Leaderboard,Spectrum,N):
    if len(Leaderboard) > N:
        results = []
        for Peptide in Leaderboard:
            try:
                Peptide_experimental_spectrum = generate_subspectrums(Peptide)
            except:
                Peptide = Peptide[0]+[Peptide[1]]
                Peptide_experimental_spectrum = generate_subspectrums(Peptide)
            results.append((Peptide,score(Spectrum,Peptide_experimental_spectrum)))
        tie = heapq.nlargest(N,results,key=lambda x: x[1])[-1][1]
        res = list(ifilter(lambda x: x[1]>=tie,results))
        return list(izip(*res))[0]
    else:
        return Leaderboard
            
def LeaderboardCyclopeptideSequencing(M,N,Spectrum):
    Leaderboard = [0]
    LeaderPeptide = []
    table = select_spectrum(Spectrum,M)
    while Leaderboard != []:
        Leaderboard = [list(pt) for pt in product(Leaderboard,table)]
        for Peptide in Leaderboard:
            try:
                Peptide_experimental_spectrum = generate_subspectrums(Peptide)
            except:
                Leaderboard = [Peptide[0]+[Peptide[1]] if x == Peptide else x for x in Leaderboard]
                Peptide = Peptide[0]+[Peptide[1]]
                Peptide_experimental_spectrum = generate_subspectrums(Peptide)
            if max(Peptide_experimental_spectrum) == max(Spectrum):
                LeaderPeptide_experimental_spectrum = generate_subspectrums(LeaderPeptide)
                if score(Spectrum,Peptide_experimental_spectrum) > score(Spectrum,LeaderPeptide_experimental_spectrum):
                    LeaderPeptide = Peptide
            elif max(Peptide_experimental_spectrum) > max(Spectrum):
                Leaderboard.remove(Peptide)
        Leaderboard = cut(Leaderboard,Spectrum,N)
    return LeaderPeptide

def result(filename):
    M,N,Spectrum = read_file(filename)
    results = LeaderboardCyclopeptideSequencing(M,N,Spectrum)
    return results[1:]

if __name__ == "__main__":

    start = timeit.default_timer()
    results = result(sys.argv[-1])
    print '-'.join(map(str,results))
    print ''
    stop = timeit.default_timer()
    print stop - start
