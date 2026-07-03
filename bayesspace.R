library(Seurat)
library(dplyr)
library(openxlsx)
library(BayesSpace)
sce <- readVisium('/BGFS1/projectdata/project_anzhen_liyulin/Cellranger/SPATIAL/IR72h/outs/')
set.seed(1000)
sce <- spatialPreprocess(sce, platform="Visium", 
                         n.PCs=15, n.HVGs=2000)
sce <- qTune(sce, qs=seq(2,20), platform="Visium", d=15)
qPlot(sce)

set.seed(1000)
sce <- spatialCluster(sce, q=11, platform="Visium", d=15,
                      init.method="mclust", model="t", gamma=3,
                      save.chain=TRUE)

write.table(colData(sce), 'bs_meta.txt', quote = F, sep = '\t')


seurat <- Seurat::CreateSeuratObject(counts=logcounts(sce),
                                     assay='Spatial',
                                     meta.data=as.data.frame(colData(sce)))
seurat <- Seurat::SetIdent(seurat, value = "spatial.cluster")
saveRDS(sce, 'IR72h_bayesspace.rds')
