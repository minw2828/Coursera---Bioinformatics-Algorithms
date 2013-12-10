#!/usr/bin/python

###############################################################################
#
# Author:
#
# Sanyk28 (san-heng-yi-shu@163.com)
#
# Date created:
#
# 4 Dec 2013
#
# CODE CHALLENGE: Solve the Eulerian Cycle Problem.
#      Input: The adjacency list of an Eulerian directed graph.
#      Output: An Eulerian cycle in this graph.
#
# Sample Input:
#      0 -> 3
#      1 -> 0
#      2 -> 1,6
#      3 -> 2
#      4 -> 2
#      5 -> 4
#      6 -> 5,8
#      7 -> 9
#      8 -> 7
#      9 -> 6
#
# Sample Output:
#      6->8->7->9->6->5->4->2->1->0->3->2->6
#
############################################################################### 

import sys
import timeit
import re
import heapq
import random
import numpy as np
from scipy import stats
from itertools import combinations,product,izip,ifilter,chain
from collections import Counter,defaultdict

def read_file(input_file):
    '''
    >>> data = read_file('test.eulerian_cycle.txt')
    >>> data
    [('0', '3'), ('1', '0'), ('2', '1'), ('2', '6'), ('3', '2'), ('4', '2'), ('5', '4'), ('6', '5'), ('6', '8'), ('7', '9'), ('8', '7'), ('9', '6')]
    >>> data = read_file('test.eulerian_cycle.3.txt')
    >>> data = read_file('test.eulerian_cycle.extra.txt')
    '''
    f = open(input_file)
    data = [item.strip() for item in f.readlines()]
    f.close()
    result = []
    for line in data:
        dp = data_process(line)
        if isinstance(dp,tuple):
            result.append(dp)
        elif isinstance(dp,list):
            result += dp
    return result

def data_process(line):
    line = [item.strip() for item in line.split('->')]
    if ',' in line[1]:
        for item in line[1].split(','):
            line.append((line[0],item))
        line.pop(0)
        line.pop(0)
    else:
        line = tuple(line)
    return line

def form_cycle(data):
    rgk,rgv = random.choice(data)
    cycle = [rgk]
    while len(data) != 0:
        try:
            cycle.append(rgv)
            data.remove((rgk,rgv))
            rgk = rgv
            rgv = [item2 for (item1,item2) in data if rgv == item1][0]
        except:
            break
    return (cycle,data)

def form_unCycle(data,rgk):
    rgv = [item2 for (item1,item2) in data if rgk == item1][0]
    cycle = [rgk]
    while len(data) != 0:
        try:
            cycle.append(rgv)
            data.remove((rgk,rgv))
            rgk = rgv
            rgv = [item2 for (item1,item2) in data if rgv == item1][0]
        except:
            break
    return (cycle,data)

def fuse(Cycle,new_Cycle):
    fuse_index = Cycle.index(new_Cycle[0])
    return Cycle[:fuse_index]+new_Cycle+Cycle[fuse_index+1:]

def eulerian_cycle(data):
    Cycle,unCycle = form_cycle(data)
    while len(unCycle) != 0:
        potential = [item for item in Cycle if item in chain(*unCycle)]
        newStart = random.choice(potential)
        new_Cycle,unCycle = form_unCycle(unCycle,newStart)
        Cycle = fuse(Cycle,new_Cycle)
    return Cycle

def result(filename):
    data = read_file(filename)
    results = eulerian_cycle(data)
    return results

if __name__ == "__main__":

    start = timeit.default_timer()
    results = result(sys.argv[-1])
    fw = open('output.'+sys.argv[-1][:-4]+'.txt','w')
    fw.write('->'.join(results))
    fw.close()
    stop = timeit.default_timer()
    print stop - start
