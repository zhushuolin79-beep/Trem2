#! /usr/bin/env python

from collections import defaultdict
from collections import Counter

clust_lr = defaultdict(list)
with open('spatial_lr_sig_info.txt', 'r') as handle:
    for line in handle.readlines():
        cell_id, cluster, lr_pair, pval = line.strip().split(' ')
        clust_lr[cluster].append(lr_pair)


with open ('../Trem2_MYO_EC_lr_res.txt', 'r') as handle:
    title = handle.readline().strip()
    print ('%s\tIn_Sp_cluster' %(title))
    for line in handle.readlines():
        sender, ligand, receiver, receptor, lr_pair, weight, mean, pval = line.strip().split('\t')
        lrs_sp = clust_lr['spa_gcn_5']
        lr_sc = ligand.capitalize() + '_' + receptor.capitalize()
        if lr_sc in lrs_sp:
            print ('%s\t%s' %(line.strip(), 'Yes'))
        else:
            print ('%s\t%s' %(line.strip(), 'No'))
