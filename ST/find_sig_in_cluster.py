#! /usr/bin/env python

cell_ids = []
cell_clust = {}
with open('IR72_meta.txt', 'r') as handle:
    title = handle.readline()
    for line in handle.readlines():
        cell_id, *_, refined_pred = line.strip().split('\t')
        cell_ids.append(cell_id)
        cell_clust[cell_id] = 'spa_gcn_' + refined_pred


lrs = []
with open('ir72_spatial_lr_sum.txt', 'r') as handle:
    title = handle.readline()
    for line in handle.readlines():
        lr_pair, *_ = line.strip().split('\t')
        lrs.append(lr_pair)


with open('ir72_lr_pval.txt', 'r') as handle:
    title = handle.readline().strip()
    # print ('\t%s' %('\t'.join(lrs)))
    lines = handle.readlines()
    for i in range(0, len(cell_ids)):
        line = lines[i]
        cell_id = cell_ids[i]
        index, *padjs = line.strip().split('\t')
        for n in range(0, len(padjs)):
            lr_pair = lrs[n]
            padj = padjs[n]
            if float(padj) <= 0.05:
                print (cell_id, cell_clust[cell_id], lr_pair, padj)
        # print ('%s\t%s' %(cell_id, '\t'.join(scores)))
