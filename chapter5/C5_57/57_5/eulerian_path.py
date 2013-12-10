#!/usr/bin/python

###############################################################################
#
# Author:
#
# Sanyk28 (san-heng-yi-shu@163.com)
#
# Date created:
#
# 5 Dec 2013
#
# CODE CHALLENGE: Solve the Eulerian Path Problem.
#      Input: The adjacency list of a directed graph that has an Eulerian path.
#      Output: An Eulerian path in this graph.
#
# Sample Input:
#      0 -> 2
#      1 -> 3
#      2 -> 1
#      3 -> 0,4
#      6 -> 3,7
#      7 -> 8
#      8 -> 9
#      9 -> 6
#
# Sample Output:
#      6->7->8->9->6->3->0->2->1->3->4
#
############################################################################### 

import sys
import timeit
import re
import heapq
import random
import numpy as np
from itertools import combinations,product,izip,ifilter,chain
from collections import Counter,defaultdict

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

def read_file(input_file):
    '''
    >>> data = read_file('test.eulerian_path.txt')
    >>> data
    [('0', '2'), ('1', '3'), ('2', '1'), ('3', '0'), ('3', '4'), ('6', '3'), ('6', '7'), ('7', '8'), ('8', '9'), ('9', '6')]
    >>> data = read_file('test.eulerian_path.extra.txt')
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

def find_ender(data):
    item0 = [item[0] for item in data]
    item1 = [item[1] for item in data]
    c0 = Counter(item0)
    c1 = Counter(item1)
    result = [(item,c0[item],c1[item]) for item in c1 if c0[item] != c1[item]]
    return [item[0] for item in result if item[1] < item[2]][0]

def find_path(data):
    initial_ender = find_ender(data)
    ender = find_ender(data)
    path = [ender]
    while len(data) != 0:
        try:
            starter = [item1 for (item1,item2) in data if item2 == ender][0]
            data.remove((starter,ender))
            path.append(starter)
            ender = starter
        except:
            break
    return (path[::-1],data,initial_ender)

def form_cycle(unpath,initial_ender):
    rgk,rgv = [(item1,item2) for (item1,item2) in unpath if item2 == initial_ender][0]
    cycle = [rgv]
    while len(unpath) != 0:
        try:
            cycle.append(rgk)
            unpath.remove((rgk,rgv))
            rgk,rgv = [(item1,item2) for (item1,item2) in unpath if item2 == rgk][0]
        except:
            break
    return (cycle[::-1],unpath)

def fuse(path,cycle):
    fuse_index = path.index(cycle[0])
    return path[:fuse_index]+cycle+path[fuse_index+1:]

def eulerian_path(data):
    path,unpath,initial_ender = find_path(data)
    cycle,unpath = form_cycle(unpath,initial_ender)
    path = fuse(path,cycle)
    while len(unpath) != 0:
        print 'len(unpath): '+str(len(unpath))
        potential = [item for item in path if item in chain(*unpath)]
        newStart = random.choice(potential)
        cycle,unpath = form_unCycle(unpath,newStart)
        path = fuse(path,cycle)   
    return path

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

def result(filename):
    data = read_file(filename)
    results = eulerian_path(data)
    return results

if __name__ == "__main__":

    start = timeit.default_timer()
    results = result(sys.argv[-1])
    fw = open('output.'+sys.argv[-1][:-4]+'.txt','w')
    fw.write('->'.join(results))
    fw.close()
    stop = timeit.default_timer()
    print stop - start
