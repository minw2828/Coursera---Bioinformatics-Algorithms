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
# CODE CHALLENGE: Implement LEADERBOARDCYCLOPEPTIDESEQUENCING.
#
# Input: Integer N and a collection of integers Spectrum.
# 
# Output: LeaderPeptide after running LEADERBOARDCYCLOPEPTIDESEQUENCING(Spectrum, N).
#
# Sample Input:
#     10
#     0 71 113 129 147 200 218 260 313 331 347 389 460
#
# Sample Output:
#     113-147-71-129
#
############################################################################### 

import sys
import timeit
import heapq
from itertools import product,izip,ifilter
from collections import Counter

table = [57, 71, 87, 97, 99, 101, 103, 113, 113, 114, 115, 128, 128, 129, 131, 137, 147, 156, 163, 186]

def read_file(input_file):
    f = open(input_file)
    N,Spectrum = [item.strip() for item in f.readlines()]
    f.close()
    return (int(N),map(int,Spectrum.split(' ')))

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
            
def LeaderboardCyclopeptideSequencing(Spectrum,N):
    Leaderboard = [0]
    LeaderPeptide = []
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
    N,Spectrum = read_file(filename)
    results = LeaderboardCyclopeptideSequencing(Spectrum,N)
    return results[1:]

if __name__ == "__main__":

    start = timeit.default_timer()
    results = result(sys.argv[-1])
    print '-'.join(map(str,results))
    print ''
    stop = timeit.default_timer()
    print stop - start
