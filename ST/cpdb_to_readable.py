#! /usr/bin/env python

from collections import defaultdict


lr_celltype_means = defaultdict(dict)
with open('out/means.txt', 'r') as handle:
    title = handle.readline().strip()
    infos = title.split('\t')
    interacting_celltypes = infos[11:]
    for line in handle.readlines():
        temp = line.strip().split('\t')
        lr_pair = temp[1]
        mean_values = temp[11:]
        for i in range(0, len(interacting_celltypes)):
            lr_celltype_means[lr_pair][interacting_celltypes[i]] = mean_values[i]

lr_celltype_pvals = defaultdict(dict)
with open('out/pvalues.txt', 'r') as handle:
    title = handle.readline().strip()
    infos = title.split('\t')
    interacting_celltypes = infos[11:]
    for line in handle.readlines():
        temp = line.strip().split('\t')
        lr_pair = temp[1]
        pvals = temp[11:]
        for i in range(0, len(interacting_celltypes)):
            lr_celltype_pvals[lr_pair][interacting_celltypes[i]] = pvals[i]

print ('LR_pair\tInter\tmeans\tpval')
for lr in lr_celltype_means:
    for inter in lr_celltype_means[lr]:
        means = lr_celltype_means[lr][inter]
        pvals = lr_celltype_pvals[lr][inter]
        if float(pvals) < 0.05:
            print ('%s\t%s\t%s\t%s' %(lr, inter, means, pvals))
