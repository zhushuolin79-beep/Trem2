#! /usr/bin/env python

from collections import Counter
from collections import defaultdict

lr_clust = defaultdict(list)
with open('spatial_lr_sig_info.txt', 'r') as handle:
    for line in handle.readlines():
        cell, cluster, lr_pair, pval = line.strip().split(' ')
        if cluster in ['spa_gcn_0', 'spa_gcn_2', 'spa_gcn_4', 'spa_gcn_5']:
            lr_clust[lr_pair].append(cluster)


for pair in lr_clust:
    nums = []
    for cluster in sorted(Counter(lr_clust[pair])):
        num = Counter(lr_clust[pair])[cluster]
        nums.append(cluster + ' : ' + str(num))
    print ('%s\t%s' %(pair, '\t'.join(nums)))
    