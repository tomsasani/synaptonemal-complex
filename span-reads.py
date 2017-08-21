import sys
import toolshed as ts
from collections import defaultdict

tsv_dict = ts.reader(sys.argv[1], header='ordered')

hq_read_dict = defaultdict(list)
hq_read_meth = defaultdict(list)
min_r = 5.

for site in tsv_dict:
    if float(site['log_lik_ratio']) > min_r:
        hq_read_dict[site['read_name']].extend([int(site['start']), int(site['end'])])
        hq_read_meth[site['read_name']].append((site['chromosome'], site['start'], site['end']))

for read in hq_read_dict:
    if max(hq_read_dict[read]) - min(hq_read_dict[read]) > 30000:
        meth_list = hq_read_meth[read]
        sorted_list = sorted(meth_list, key=lambda tup: tup[1])
        if len(set([x[0] for x in sorted_list])) > 1:
            continue
        #print '\t'.join([sorted_list[0][0], sorted_list[0][1], sorted_list[-1][2]])
        for site in sorted_list:
            """
            if sorted_list.index(site) == 0:
                print '\t'.join([site[0], site[1], site[2], '-100'])
            elif sorted_list.index(site) == len(sorted_list) - 1:
                print '\t'.join([site[0], site[1], site[2], '100'])
            else: 
                print '\t'.join([site[0], site[1], site[2], '1'])
            """
            print '\t'.join([site[0], site[1], site[2]])




