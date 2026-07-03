import os,csv,re
import pandas as pd
import numpy as np
import scanpy as sc
import math
import SpaGCN as spg
from scipy.sparse import issparse
import random, torch
import warnings
warnings.filterwarnings("ignore")
import matplotlib.colors as clr
import matplotlib.pyplot as plt
import SpaGCN as spg
import cv2

adata = sc.read_visium('/BGFS1/projectdata/project_anzhen_liyulin/Cellranger/SPATIAL/IR72h/outs')
x_array = adata.obs['array_row'].tolist()
y_array = adata.obs['array_col'].tolist()
x_pixel = [item[1] for item in adata.obsm['spatial']]
y_pixel = [item[0] for item in adata.obsm['spatial']]
img = cv2.imread('IR72h.V10J29-061.A1.tif')

s=1
b=49
adj=spg.calculate_adj_matrix(x=x_pixel,y=y_pixel, x_pixel=x_pixel, y_pixel=y_pixel, image= img, beta=b, alpha=s, histology=False)

adata.var_names_make_unique()
spg.prefilter_genes(adata,min_cells=3)
spg.prefilter_specialgenes(adata)
sc.pp.normalize_total(adata, target_sum=10000)
sc.pp.log1p(adata)

p=0.5 
l=spg.search_l(p, adj, start=0.01, end=1000, tol=0.01, max_run=100)

n_clusters=7
r_seed=t_seed=n_seed=100
res=spg.search_res(adata, adj, l, n_clusters, start=0.7, step=0.1, tol=5e-3, lr=0.05, max_epochs=20, r_seed=r_seed, t_seed=t_seed, n_seed=n_seed)

clf=spg.SpaGCN()
clf.set_l(l)
#Set seed
random.seed(r_seed)
torch.manual_seed(t_seed)
np.random.seed(n_seed)
#Run
clf.train(adata,adj,init_spa=True,init="louvain",res=res, tol=5e-3, lr=0.05, max_epochs=200)
y_pred, prob=clf.predict()
adata.obs["pred"]= y_pred
adata.obs["pred"]=adata.obs["pred"].astype('category')
adj_2d=spg.calculate_adj_matrix(x=x_array,y=y_array, histology=False)
refined_pred=spg.refine(sample_id=adata.obs.index.tolist(), pred=adata.obs["pred"].tolist(), dis=adj_2d, shape="hexagon")
adata.obs["refined_pred"]=refined_pred
adata.obs["refined_pred"]=adata.obs["refined_pred"].astype('category')
adata.write_h5ad("./new_results.h5ad")

plt.rcParams["figure.figsize"] = (8, 8)
for i in groups:
    sc.pl.spatial(adata, img_key='hires', color='refined_pred', groups = [int(i)], size=1.5)
