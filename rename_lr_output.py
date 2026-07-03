#! /usr/bin/env python

from math import log10

cell_ids = []
with open('IR72_meta.txt', 'r') as handle:
    title = handle.readline()
    for line in handle.readlines():
        cell_id, *_ = line.strip().split('\t')
        cell_ids.append(cell_id)


lrs = []
with open('ir72_spatial_lr_sum.txt', 'r') as handle:
    title = handle.readline()
    for line in handle.readlines():
        lr_pair, *_ = line.strip().split('\t')
        lrs.append(lr_pair)

with open('ir72_lr_pval.txt', 'r') as handle:
    title = handle.readline().strip()
    print ('\t%s' %('\t'.join(lrs)))
    lines = handle.readlines()
    for i in range(0, len(cell_ids)):
        line = lines[i]
        cell_id = cell_ids[i]
        index, *scores = line.strip().split('\t')
        new_value = [str(-log10(float(x))) for x in scores ]
        print ('%s\t%s' %(cell_id, '\t'.join(new_value)))