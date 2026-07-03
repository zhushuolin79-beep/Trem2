#! /usr/bin/env python

from collections import defaultdict
from collections import Counter

ligs = []
recs = []
lr_to_pval = defaultdict(dict)
lr_to_sp = defaultdict(dict)


with open('sc_lr_spa_0_info.txt', 'r') as handle:
    title = handle.readline().strip()
    # print (title)
    for line in handle.readlines():
        sender, ligand, receiver, receptor, lr_pair, weight, mean, pval, sp = line.strip().split('\t')
        ligs.append(ligand)
        recs.append(receptor)
        lr_to_sp[ligand][receptor] = sp
        lr_to_pval[ligand][receptor] = pval

print ('%s\t%s' %('', '\t'.join(Counter(recs))))
for ligand in Counter(ligs):
    rec_vals = []
    for receptor in Counter(recs):
        if receptor in lr_to_pval[ligand] and lr_to_sp[ligand][receptor] == 'Yes':
            rec_vals.append('1.001')
        elif receptor in lr_to_pval[ligand] and lr_to_sp[ligand][receptor] == 'No':
            rec_vals.append('1')
        else:
            rec_vals.append('0')
    print ('%s\t%s' %(ligand, '\t'.join(rec_vals)))