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

p = argparse.ArgumentParser()
p.add_argument('--f', nargs='*')
p.add_argument('-chrom', default='all')
p.add_argument('-w', type=int, default=500)
args = p.parse_args()

a, b = defaultdict(list), defaultdict(list)

def add_to_dict(f, d):
    with open(f) as tsvfile:
        for line in csv.reader(tsvfile, delimiter='\t'):
            start, end = int(line[1]), int(line[2])
            multiplier = int(float(line[3])) 
            d[line[0]].append((start, end, multiplier))
    return d

def plot(d, w=500, chrom_specified='all', color='red'):


    for chrom in d:
        xy_vals = defaultdict(list)
        if chrom_specified != 'all' and chrom_specified != chrom:
            continue
        max_pos = max([x[1] for x in d[chrom]])
        max_mult = max([x[2] for x in d[chrom]])
        sliding_window = [(x, x + w) for x in range(1, max_pos) if x % w == 0 or x == 0]
        values = d[chrom]
        for val in values:
            window = [w for w in sliding_window if val[0] in range(w[0], w[1])]
            multiplier = val[2]
            if len(window) > 0:
                xy_vals[window[0]].append(multiplier)

        if chrom_specified != 'all':
            for xy in xy_vals:
                plt.plot([xy[0], xy[1]], [np.mean(xy_vals[xy]), np.mean(xy_vals[xy])], c=color)
        else:
            f, axes = plt.subplots(4, 4, figsize=(15, 15))
            
            lst = [0, 1, 2, 3]

            p = itertools.permutations(lst, 2)
            p = [list(x) for x in p]
            p.append([0,0])
            p.append([1,1])
            p.append([2,2])
            p.append([3,3])

            p = sorted(p, key=lambda x: (x[0], x[1]))
            idx = 0
            for xy in xy_vals:
                plt.plot([xy[0], xy[1]], [np.mean(xy_vals[xy]), np.mean(xy_vals[xy])], ax=axes[p[idx][0], p[idx][1]])
            sns.despine(top=True, left=True)
            idx += 1
        plt.savefig('out.eps')

colors = {0:'red', 1:'blue'}
count = 0
for f in args.f:
    di = defaultdict(list)
    d = add_to_dict(f, di)
    plot(d, w=args.w, chrom_specified=args.chrom, color=colors[count])
    count += 1

"""
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
"""
