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

#rc('text', usetex=True)

sns.set(style="white", palette="muted", color_codes=True)

p = argparse.ArgumentParser()
p.add_argument('--f', nargs='*')
p.add_argument('-chrom', default='all')
args = p.parse_args()

def add_to_dict(f):

    d = defaultdict(list)

    ftype = None
    if f.split('.')[-1] == 'bedgraph':
        ftype = 'bedgraph'
        with open(f) as tsvfile:
            for line in csv.reader(tsvfile, delimiter='\t'):
                start, end = int(line[1]), int(line[2])
                multiplier = int(float(line[3]))
                if multiplier < 0:
                    continue
                for i in range(multiplier):
                    d[line[0]].append(start)

    elif f.split('.')[-1] == 'bed':
        ftype = 'bed'
        with open(f) as tsvfile:
            for line in csv.reader(tsvfile, delimiter='\t'):
                start, end = int(line[1]), int(line[2])
                d[line[0]].extend((range(start, end + 1)))
    return d, ftype 

def plot(d, axes, chrom_specified='all', color='red', sample=None, ftype=None):
    
    idx = 0
    if chrom_specified == 'all':
        lst = [0, 1, 2, 3]

        p = itertools.permutations(lst, 2)
        p = [list(x) for x in p]
        p.append([0,0])
        p.append([1,1])
        p.append([2,2])
        p.append([3,3])

        p = sorted(p, key=lambda x: (x[0], x[1]))

    for chrom in d:
        if chrom == 'chrM':
            continue
        if chrom_specified != 'all' and chrom_specified != chrom:
            continue
        values = d[chrom]
        if chrom_specified != 'all':
            if ftype == 'bed':
                sns.rugplot(np.array(values), color=color)
            else:
                sns.distplot(np.array(values), color=color, rug=False, hist=False,
                                                                   hist_kws={"histtype":"step",
                                                                             "lw":1,
                                                                             "color":color},
                                                                   kde_kws={"color":color,
                                                                            "lw":2,
                                                                            "bw":0.1},
                                                                   label=sample)
        else:
            if ftype == 'bed':
                sns.rugplot(np.array(values), color=color)
            else:
                sns.distplot(np.array(values), color=color, ax=axes[p[idx][0], p[idx][1]], 
                                                                   rug=False, hist=False,
                                                                   hist_kws={"histtype":"step",
                                                                             "lw":1,
                                                                             "color":color},
                                                                   kde_kws={"color":color,
                                                                            "lw":2,
                                                                            "bw":0.1},
                                                                   label=sample)
            sns.despine(top=True, left=True)
            idx += 1

    keys = d.keys()
    keys = [x for x in keys if x != 'chrM']
    keys = sorted(keys, key=lambda x: int(x.split('chr')[1]))

    if args.chrom == 'all':
        count = 0
        for ax in axes:
            for a in ax:
                a.get_yaxis().set_visible(False)
                a.get_xaxis().set_visible(False)
            #    a.legend().set_visible(False)
                a.set_title(r"\textbf{%s}" % keys[count])
                count += 1
    plt.savefig('out.eps')

colors = sns.color_palette("hls", len(args.f))
count = 0
if args.chrom == 'all':
    fig, axes = plt.subplots(4, 4, figsize=(15, 15))
else:
    fig, axes = None, None
for f in args.f:
    d, ftype = add_to_dict(f)
    plot(d, axes, chrom_specified=args.chrom, color=colors[count], sample=f, ftype=ftype)
    count += 1
