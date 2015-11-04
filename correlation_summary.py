#!/usr/bin/python
import sys
import csv
import gzip
import numpy
from scipy.stats.stats import pearsonr

def num(s):
    try:
        return int(s)
    except ValueError:
        return float(s)

def LoadData(fnamePairs, corrs):
    M = -1
    with open(fnamePairs, 'r') as f:
        reader=csv.reader(f, delimiter='\t')
        for line in reader:
            w = True
            l = []
            if M == -1:
                M = len(line[1:])
            elif M != len(line[1:]):
                print "Unexpected input: non-equal number of entries in lines."
                sys.exit(1)
            for val in line[1:]:
                if val == "nan":
                    w = False
                l.append(num(val))
            if w:
                corrs[line[0]] = l
    return M

def CreateReport(corrs, id):
    report = "shift = " + `id`
    l = []
    for key in corrs:
        l.append(corrs[key][id])
    report += ", average correlation = " +  `sum(l)/len(l)`
    report += ", sd = " + `numpy.std(l)`
    print report
    

try:
    fname = sys.argv[1]
except IndexError:
    print "Please specify the input file"
    sys.exit(1)

#fname = "0897557-mincut-nopred.tsv.gz"
#fnamePairs = "a1c3d1e-naive-nopred.tsv.gz"
corrs = dict()
M = LoadData(fname, corrs)

averCor = [0]*M
for key in corrs:
    for i in range(M-1):
        averCor[i] += corrs[key][1+i]
for i in range(M):
    averCor[i] = averCor[i] / len(corrs)
#corrs['average'] = averCor
#print averCor
maxCor = max(averCor)
ind = [i for i, j in enumerate(averCor) if j == maxCor]
if ind[0] != 0:
    CreateReport(corrs, 0)
for i in ind:
    CreateReport(corrs, i)