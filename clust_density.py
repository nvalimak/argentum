#!/usr/bin/python
import sys
import csv
import gzip
from scipy.stats.stats import pearsonr
import numpy
import math

def LoadData(fnamePairs, fPopName):
    with open(fPopName, 'r') as f:
        reader=csv.reader(f, delimiter='\t')
        hap2Pop = dict()
        for line in reader:
            hap2Pop[int(line[0])] = int(line[1])
    M = max(hap2Pop.keys())
    print M
    p1 = min(hap2Pop.values())
    p2 = max(hap2Pop.values())
    clustDen = dict()
    for i in range(p1, p2+1):
        for j in range(i, p2+1):
            key = `i` + '-' + `j`
            clustDen[key] = [0]*(M-1)
    
    with gzip.open(fnamePairs, 'r') as f:
        reader=csv.reader(f, delimiter='\t')
        for line in reader:
            p1 = hap2Pop[int(line[1])]
            p2 = hap2Pop[int(line[2])]
            key = `min(p1,p2)` + '-' + `max(p1,p2)`
            clustDen[key][int(line[4])-2] += 1
            clustDen[key][int(line[3])-2] -= 1
#            scrm.append( num(line[3]) ) 
#            argentum.append( num(line[4]) )
    PrintData(clustDen)



def PrintData(clustDen):
    for key in clustDen:
        line = "[" + key + "],[" + `numpy.std(clustDen[key])` + "]"
        for v in clustDen[key]:
            line = line + "," + `v`
        print line



try:
    fname = sys.argv[1]
except IndexError:
    print "Please specify the input file"
    sys.exit(1)

try:
    fPopName = sys.argv[2]
except IndexError:
    print "Please specify the population input file"
    sys.exit(1)

LoadData(fname, fPopName)