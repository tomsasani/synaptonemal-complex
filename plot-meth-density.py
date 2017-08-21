import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rc
import seaborn as sns
import numpy as np
import scipy
import itertools
import argparse
from collections import defaultdict
import sys
import csv

rc('text', usetex=True)

sns.set(style="white", palette="muted", color_codes=True)

a, b, c = defaultdict(list), defaultdict(list), defaultdict(list)

def add_to_dict(f, d):
    if f.split('.')[-1] == 'bedgraph':
        with open(f) as tsvfile:
            for line in csv.reader(tsvfile, delimiter='\t'):
                start, end = int(line[1]), int(line[2])
                multiplier = int(float(line[3])) 
                if multiplier < 0:
                    continue
                for i in range(multiplier):
                    d[line[0]].append(start)
    elif f.split('.')[-1] == 'bed':
        with open(f) as tsvfile:
            for line in csv.reader(tsvfile, delimiter='\t'):
                start, end = int(line[1]), int(line[2])
                d[line[0]].extend(range(start, end + 1))

    return d

red1, ont = add_to_dict(sys.argv[2], a), add_to_dict(sys.argv[1], b)

lst = [0, 1, 2, 3]

p = itertools.permutations(lst, 2)
p = [list(x) for x in p]
p.append([0,0])
p.append([1,1])
p.append([2,2])
p.append([3,3])

p = sorted(p, key=lambda x: (x[0], x[1]))
pvals = defaultdict()
keys = red1.keys()
keys = sorted(keys, key=lambda x: int(x.split('chr')[1]))

f, axes = plt.subplots(4, 4, figsize=(15, 15))
idx = 0
for chrom in keys:
    if len(ont[chrom]) > 0:
        sns.distplot(np.array(ont[chrom]), hist=True, ax=axes[p[idx][0], p[idx][1]], color='blue', label='ont', norm_hist=False)
    if len(red1[chrom]) > 0:
        sns.distplot(np.array(red1[chrom]), hist=True, ax=axes[p[idx][0], p[idx][1]], color='red', label='elife', norm_hist=False)
    if len(ont[chrom]) > 0 and len(red1[chrom]) > 0:
        _, pval = scipy.stats.ks_2samp(np.array(ont[chrom]), np.array(red1[chrom]))
        pvals[idx] = pval
    sns.despine(top=True, left=True)
    idx += 1
count = 0
for ax in axes:
    for a in ax:
        a.get_yaxis().set_visible(False)
        a.get_xaxis().set_visible(False)
        a.legend().set_visible(False)
        if pvals[count] < 0.05:
            a.set_title(r"\textbf{%s}" % keys[count])
        count += 1
plt.savefig('out.eps')
