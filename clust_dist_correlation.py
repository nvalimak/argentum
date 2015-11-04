#!/usr/bin/python
import sys
import csv
import gzip
from scipy.stats.stats import pearsonr


def num(s):
    try:
        return int(s)
    except ValueError:
        return float(s)

def CorrelationRange(r, scrm, argentum, cor12):
	for shift in range(r):
		cor12.append ( Corr(shift, scrm, argentum) )


def Corr(shift, scrm, argentum):
    l = len(scrm)
    if shift >= 0:
		b1 = shift
		b2 = 0
		e1 = l
		e2 = l - shift
    else:
        b1 = 0
        b2 = 0 - shift
        e1 = l + shift
        e2 = l
#    scrm = [row[1] for row in clusterDist]
#    argentum = [row[2] for row in clusterDist]
    return pearsonr(scrm[b1:e1], argentum[b2:e2])[0]









def LoadData(fnamePairs, r):
    corrs = dict()
    with gzip.open(fnamePairs, 'r') as f:
        reader=csv.reader(f, delimiter='\t')#, fieldnames=['site', 'leaf_1', 'leaf_2', 'scrv', 'argentum'])
        scrm = []
        argentum = []
        l1 = 0
        l2 = 0
        for line in reader:
            if l1 == l2:
                l1 = line[1]
                l2 = line[2]
            elif l1 != line[1] or l2 != line[2]:
                cor12 = []
                CorrelationRange(r, scrm, argentum, cor12)
                entry = l1 + "," + l2
                corrs[entry] = cor12
                scrm = []
                argentum = []
                l1 = line[1]
                l2 = line[2]
            scrm.append( num(line[3]) ) 
            argentum.append( num(line[4]) )
            
        cor12 = []
        CorrelationRange(r, scrm, argentum, cor12)
        entry = l1 + "," + l2
        corrs[entry] = cor12
    for key in corrs:
        line = key
        for v in corrs[key]:
            line = line + "\t" + `v`
        print line

try:
    fname = sys.argv[1]
except IndexError:
    print "Please specify the input file"
    sys.exit(1)

try:
    r = num(sys.argv[2]) #Second parameter is optional. It corresponds to the range of shift [0, r)
except IndexError:
    r = 0 

#fname = "0897557-mincut-nopred.tsv.gz"
#fname = "a1c3d1e-naive-nopred.tsv.gz"
LoadData(fname, r)

#zcat < a1c3d1e-naive-nopred.tsv.gz | sort -n -k1,1n -k2,2n | more > a1c3d1e-naive-nopred-sorted.tsv