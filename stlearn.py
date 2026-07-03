import stlearn as st
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
from collections import Counter
#分析原始数据
data_dir = 'E:/ST/analysis/STdata/data/IR0h/outs/'
ir72 = st.Read10X(data_dir)
st.add.image(adata=ir72,
             imgpath=data_dir+"spatial/tissue_hires_image.png",
             library_id="IR0h", visium=True)
ir72.var_names_make_unique()
st.pp.filter_genes(ir72, min_cells=3)
st.pp.normalize_total(ir72)

#添加zone列

#Running the Ligand-Receptor Analysis
lrs = st.tl.cci.load_lrs(['connectomeDB2020_lit'], species='mouse')


cust_lr = []
with open('/BGFS1/home/fengk/project/anzhen/IR72_nichenet/5_6_newmarker/lrs_for_spatial.txt', 'r') as handle:
    title = handle.readline().strip()
    for line in handle.readlines():
        lr = line.strip()
#         lig, lr_pair, weight = line.strip().split('\t')
#         ligand, receptor = lr_pair.split('_')
#         new_lr = ligand.capitalize() + '_' + receptor.capitalize()
        cust_lr.append(lr)
cust_lr = list(set(cust_lr))

lrs = cust_lr
st.tl.cci.run(ir72, lrs,
                  min_spots = 20, #Filter out any LR pairs with no scores for less than min_spots
                  distance=None, # None defaults to spot+immediate neighbours; distance=0 for within-spot mode
                  n_pairs=10000, # Number of random pairs to generate; low as example, recommend ~10,000
                  n_cpus=16, # Number of CPUs for parallel. If None, detects & use all available.
                  )
                  
                  
ir72.uns['lr_summary'].to_csv('ir72_spatial_lr_sum.txt', sep = '\t')

ir72.write_h5ad("./ir72_spatial_LR.h5ad")

ir72_lr_pval = pd.DataFrame(ir72.obsm['p_vals'])
ir72_lr_padj = pd.DataFrame(ir72.obsm['p_adjs'])
ir72_lr_score = pd.DataFrame(ir72.obsm['lr_scores'])

plt.rcParams["figure.figsize"] = (10, 10)
plt.rcParams['figure.dpi'] = 100
best_lr = 'Gas6_Mertk'
stats = ['lr_scores', '-log10(p_adjs)']
st.pl.lr_result_plot(ir72, use_result='lr_scores', use_lr=best_lr, show_color_bar=True, size=50)
st.pl.lr_result_plot(ir72, use_result='-log10(p_adjs)', use_lr=best_lr, show_color_bar=True, size=50)
st.pl.lr_result_plot(ir72, use_result='lr_sig_scores', use_lr=best_lr, show_color_bar=True, size=50)


