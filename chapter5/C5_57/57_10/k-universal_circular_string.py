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
# CODE CHALLENGE: Solve the k-Universal Circular String Problem.
#      Input: An integer k.
#      Output: A k-universal circular string.
#
# Sample Input:
#      4
#
# Sample Output:
#      0000110010111101
#
############################################################################### 

import sys
import timeit
import random
from itertools import combinations,product,izip,ifilter,chain
from collections import Counter,defaultdict

def read_file(input_file):
    '''
    >>> k = read_file('test.k-universal_circular_string.4.txt')
    '''
    f = open(input_file)
    data = f.read().strip()
    f.close()
    return int(data)

def generate_binary(k):
    return [bin(i)[2:].zfill(k) for i in range(2**k)]

def form_relation(k):
    nodes = generate_binary(k-1)
    d = {}
    for item in nodes:
        add1 = item[1:]+'0'
        add2 = item[1:]+'1'
        d[item] = [add1,add2]
    return d

def form_cycle(data):
    rgk = random.choice(data.keys())
    rgv = random.choice(data[rgk])
    cycle = [rgk]
    while len(data) != 0:
        try:
            cycle.append(rgv)
            if len(data[rgk]) > 1:
                data[rgk].remove(rgv)
            else:
                del data[rgk]
            rgk = rgv
            rgv = random.choice(data[rgk])
            if rgv == cycle[0] and rgv == cycle[-1]:
                break
        except:
            break
    return (cycle,data)

def form_unCycle(data,rgk):
    choose = data[rgk]
    rgv = random.choice(choose)
    cycle = [rgk]
    while len(data) != 0:
        try:
            cycle.append(rgv)
            if len(data[rgk]) > 1:
                data[rgk].remove(rgv)
            else:
                del data[rgk]
            rgk = rgv
            rgv = random.choice(data[rgk])
            if rgv == cycle[0] and rgv == cycle[-1]:
                break
        except:
            break
    return (cycle,data)

def fuse(Cycle,new_Cycle):
    fuse_index = Cycle.index(new_Cycle[0])
    return Cycle[:fuse_index]+new_Cycle+Cycle[fuse_index+1:]

def eulerian_cycle(data,k):
    Cycle,unCycle = form_cycle(data)
    while len(unCycle) != 0:
        keys = unCycle.keys()
        potential = list(set(Cycle)&set(keys))
        newStart = random.choice(potential)
        new_Cycle,unCycle = form_unCycle(unCycle,newStart)
        Cycle = fuse(Cycle,new_Cycle)
    return Cycle

def form_string(path):
    string = ''.join([item[-1] for item in path[1:]])
    return string

def result(filename):
    k = read_file(filename)
    data = form_relation(k)
    path = eulerian_cycle(data,k)
    results = form_string(path)
    return results

if __name__ == "__main__":

    start = timeit.default_timer()
    results = result(sys.argv[-1])
    fw = open('output.'+sys.argv[-1][:-4]+'.txt','w')
    fw.write(results)
    fw.close()
    stop = timeit.default_timer()
    print stop - start
