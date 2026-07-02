rm(list = ls()) 
library(Seurat)
library(cowplot)
library(patchwork)
library(dplyr)
library(ggplot2)
library(CellChat)
library(pheatmap)
library(tidyverse)
library(tidyr)
####心脏####
setwd("D:/analysis/公共数据/")
sham <- Read10X('GSE163465_心脏/数据/sham/')
sham <- CreateSeuratObject(counts = sham , 
                            project = "sham", 
                            min.cells = 30, 
                            min.features = 500)
sham[["percent.mt"]] <- PercentageFeatureSet(sham, pattern = "^mt-")
VlnPlot(sham, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)
MI3D<- Read10X('GSE163465_心脏/Day3/')
MI3D <- CreateSeuratObject(counts = MI3D, 
                            project = "MI_3D", 
                            min.cells = 30, 
                            min.features = 500)
MI3D[["percent.mt"]] <- PercentageFeatureSet(MI3D, pattern = "^mt-")
VlnPlot(MI3D, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)
MI7D <- Read10X('GSE163465_心脏/Day7/')
MI7D <- CreateSeuratObject(counts = MI7D, 
                            project = "MI_7D", 
                            min.cells = 30, 
                            min.features = 500)
MI7D[["percent.mt"]] <- PercentageFeatureSet(MI7D, pattern = "^mt-")
VlnPlot(MI7D, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)
MI14D <- Read10X('GSE163465_心脏/Day14/')
MI14D <- CreateSeuratObject(counts = MI14D, 
                            project = "MI_14D", 
                            min.cells = 30, 
                            min.features = 500)
MI14D[["percent.mt"]] <- PercentageFeatureSet(MI14D, pattern = "^mt-")
VlnPlot(MI14D, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)

sham;MI3D;MI7D;MI14D
split_seurat <- list(sham,MI3D,MI7D,MI14D)
####心脏2####
setwd("D:/analysis/公共数据/")
sham <- Read10X('GSE210159_心脏2/sham/')
sham <- CreateSeuratObject(counts = sham , 
                           project = "sham", 
                           min.cells = 30, 
                           min.features = 500)
sham[["percent.mt"]] <- PercentageFeatureSet(sham, pattern = "^mt-")
VlnPlot(sham, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)
sham <- subset(sham, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 20)

MI1D<- Read10X('GSE210159_心脏2/MI 24H/')
MI1D <- CreateSeuratObject(counts = MI1D, 
                           project = "MI_1D", 
                           min.cells = 30, 
                           min.features = 500)
MI1D[["percent.mt"]] <- PercentageFeatureSet(MI1D, pattern = "^mt-")
VlnPlot(MI1D, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)
MI1D <- subset(MI1D, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 20)

MI3D <- Read10X('GSE210159_心脏2/MI 72H/')
MI3D <- CreateSeuratObject(counts = MI3D, 
                           project = "MI_3D", 
                           min.cells = 30, 
                           min.features = 500)
MI3D[["percent.mt"]] <- PercentageFeatureSet(MI3D, pattern = "^mt-")
VlnPlot(MI3D, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)
MI3D <- subset(MI3D, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 20)


sham;MI1D;MI3D
split_seurat <- list(sham,MI1D,MI3D)
####肾脏####
rm(list = ls())
setwd("D:/analysis/公共数据/")
data <- Read10X('GSE267242_肾脏/数据/')
data <- CreateSeuratObject(counts = data , 
                           project = "kidney", 
                           min.cells = 30, 
                           min.features = 500)
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^mt-")
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)
data <- subset(data, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 30)
#拆分样本
data@meta.data$rowLeo <- rownames(data@meta.data)
data@meta.data$Sample <- str_split(data$rowLeo,"-",simplify=T)[,2]
head(data)
split_seurat <- SplitObject(data,split.by='Sample')
split_seurat

split_seurat <- split_seurat[-3]

####肾脏2####
setwd("D:/analysis/公共数据/")
data <- readRDS("GSE180420_肾脏2/data/GSE180420_EXPORT_counts.rds")
cell <- data@Dimnames[[2]]
#control --norm1,2,3,4
#long1D -- 1,5
#short1D--2,6
#long3D--IRF4,IRF7
#short3D--IRF3,IRF11
cell_norm <- grep(paste0("^", "IRF11"), cell, value = TRUE)
data_norm <- data[,cell_norm]
data12 <- CreateSeuratObject(counts = data_norm , 
                           project = "short3D",
                           min.cells = 30, 
                           min.features = 500)
data12[["percent.mt"]] <- PercentageFeatureSet(data12, pattern = "^mt-")
VlnPlot(data12, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)

data1;data2;data3;data4;data5;data6;data7;data8;data9;data10;data11;data12

split_seurat <- list(data1,data2,data3,data4,data5,data6,data7,data8,data9,data10,data11,data12)
####下肢####
rm(list = ls())
setwd("D:/analysis/公共数据/")
data <- Read10X('GSE227075_下肢/data/')
data <- CreateSeuratObject(counts = data , 
                           project = "limbs", 
                           min.cells = 30, 
                           min.features = 500)
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^mt-")
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)
data <- subset(data, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 20)
#拆分样本
data@meta.data$rowLeo <- rownames(data@meta.data)
data@meta.data[6:9] <- str_split(data$rowLeo,"_",simplify=T)
head(data)
colnames(data@meta.data)[6:9] <- c("mouse","moudle","samlpe","barcode") 
data@meta.data$sample <- paste0(data@meta.data$mouse,"_",data@meta.data$moudle,"_",data$samlpe)
data@meta.data <- data@meta.data[,-c(5:9)]
split_seurat <- SplitObject(data,split.by='sample')
split_seurat

split_seurat <- split_seurat[9:16]

####脑####
rm(list = ls())
setwd("D:/analysis/公共数据/")
data <- read.table('GSE210986_脑/data/GSM6443690_sham.counts.tsv',header=T, sep="\t")
rownames(data) <- data$gene
data <- data[,-1]
sham <- CreateSeuratObject(counts = data, 
                           project = "sham", 
                           min.cells = 30, 
                           min.features = 500)
sham[["percent.mt"]] <- PercentageFeatureSet(sham, pattern = "^mt-")
VlnPlot(sham, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)
sham <- subset(sham, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 30)

data <- read.table('GSE210986_脑/data/GSM6443691_MCAO.counts.tsv',header=T, sep="\t")
rownames(data) <- data$gene
data <- data[,-1]
MCAO <- CreateSeuratObject(counts = data , 
                           project = "MCAO", 
                           min.cells = 30, 
                           min.features = 500)
MCAO[["percent.mt"]] <- PercentageFeatureSet(MCAO, pattern = "^mt-")
VlnPlot(MCAO, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)
MCAO <- subset(MCAO, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 30)


sham;MCAO
split_seurat <- list(sham,MCAO)
####骨髓####
rm(list = ls())
setwd("D:/analysis/公共数据/")
sham <- Read10X('GSE157244_骨髓/data/D0/')
D0 <- CreateSeuratObject(counts = sham , 
                           project = "D0", 
                           min.cells = 30, 
                           min.features = 500)
D0[["percent.mt"]] <- PercentageFeatureSet(D0, pattern = "^mt-")
VlnPlot(D0, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)
D0 <- subset(D0, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 20)

sham <- Read10X('GSE157244_骨髓/data/D1/')
D1 <- CreateSeuratObject(counts = sham , 
                         project = "D1", 
                         min.cells = 30, 
                         min.features = 500)
D1[["percent.mt"]] <- PercentageFeatureSet(D1, pattern = "^mt-")
VlnPlot(D1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)
D1 <- subset(D1, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 20)

sham <- Read10X('GSE157244_骨髓/data/D2/')
D2 <- CreateSeuratObject(counts = sham , 
                         project = "D2", 
                         min.cells = 30, 
                         min.features = 500)
D2[["percent.mt"]] <- PercentageFeatureSet(D2, pattern = "^mt-")
VlnPlot(D2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)
D2 <- subset(D2, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 20)
D0;D1;D2
split_seurat <- list(D0,D1,D2)
####样本合并####
library(future)
library(harmony)
options(future.globals.maxSize = 8000 * 1024^2)
# normalize and identify variable features for each dataset independently
ifnb.list <- lapply(X = split_seurat, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})


features <- SelectIntegrationFeatures(object.list = ifnb.list)
immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features)

# this command creates an 'integrated' data assay
immune.combined <- IntegrateData(anchorset = immune.anchors)
DefaultAssay(immune.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindClusters(immune.combined, resolution = c(0.3,0.4,0.5))
Idents(immune.combined) <- "integrated_snn_res.0.3"
DimPlot(immune.combined, reduction = "umap")
saveRDS(immune.combined,"GSE180420_肾脏2/RDS/kidney_nocelltype")

####亚群定义####
rm(list = ls())
setwd("D:/analysis/公共数据/")
scRNA <- readRDS("GSE180420_肾脏2/RDS/kidney_nocelltype")
table(scRNA$orig.ident)
Idents(scRNA) <- "orig.ident"
scRNA <- subset(scRNA,ident = c("IRF10","IRF12"),invert = T)
scRNA@meta.data$sample <- factor(scRNA@meta.data$orig.ident,
                                    labels = c("long_1D","short_3D","short_1D","short_3D","long_3D","long_1D",
                                               "short_1D","long_3D","sham","sham","sham","sham"))
scRNA@meta.data$sample  <- factor(scRNA@meta.data$sample,
                                    levels = c("sham","short_1D","long_1D","short_3D","long_3D"))
Idents(scRNA) <- "sample"
DimPlot(scRNA,label = T)
saveRDS(scRNA,"GSE227075_下肢/RDS文件/climb_mac_celltype.rds")

DimPlot(scRNA,label = T,group.by = "celltype")
marker <- read.csv("D:/课题数据/ST/analysis/IR五个时间点/marker/all_marker_cosg.csv")
marker <- marker[,-1]
head(scRNA@meta.data)
Idents(scRNA) <- "integrated_snn_res.0.3"
DimPlot(scRNA,label = F,group.by = "orig.ident")
#heart
FeaturePlot(scRNA, features = c("Cd68","Adgre1","C1qc","C1qb"),label = T,ncol = 2,order = T)#2,3,5,7-Mac
FeaturePlot(scRNA, features = c("Plac8","Ly6c2","Ms4a4c","Chil3"),label = T,ncol = 2)#11-mono
FeaturePlot(scRNA, features = c("S100a8","S100a9","Lcn2","Cxcl2"),label = T,ncol = 2)#1,9-Neu
FeaturePlot(scRNA, features = c("Cd79a","Ly6d","Ebf1","Cd79b"),label = T,ncol = 2)#8-B
FeaturePlot(scRNA, features = c("Cd3g","Cd3d","Cxcr6","Trbc2"),label = T,ncol =2 )#10-T
FeaturePlot(scRNA, features = c("Gzma","Nkg7","Klra4","Ccl5"),label = T,ncol = 2)#-NK
FeaturePlot(scRNA, features = c("Fscn1","Ccl22","Ccr7","Tmem123"),label = T,ncol = 2)#
FeaturePlot(scRNA, features = marker$CF[1:4],label = T,max.cutoff = 2)#4-CF
FeaturePlot(scRNA, features = marker$EC[1:4],label = T,max.cutoff = 2)#0,12-EC
FeaturePlot(scRNA, features = marker$SMC[1:4],label = T,max.cutoff = 2)#6-SMC
FeaturePlot(scRNA, features = marker$CM[1:4],label = T,max.cutoff = 2)#14-CM
FeaturePlot(scRNA, features = marker$Mac[1:4],label = T,max.cutoff = 2)#

FeaturePlot(scRNA, features = "Pax7",label = T,max.cutoff = 2)#Muscle stem cells(MuSCs)-4,18
FeaturePlot(scRNA, features = "Top2a",label = T,max.cutoff = 2)#Cycling basal cells-8,14
FeaturePlot(scRNA, features = "Krt15",label = T,max.cutoff = 2)#Epithelial cells-9
#brain
FeaturePlot(scRNA, features = c("Agp4", "Gfap", "Mfge8", "Aldh1l1"),label = T,max.cutoff = 2)#Astrocyte-2,7
FeaturePlot(scRNA, features = "Igkc",label = T,max.cutoff = 2)#lymphocyte-15
FeaturePlot(scRNA, features = c("Flt1", "Cldn5"),label = T,max.cutoff = 2)#Endothelial cells-5
FeaturePlot(scRNA, features = "Col1a2",label = T,max.cutoff = 2)#Fibroblast-8
FeaturePlot(scRNA, features = c("Cd68","Trem2","Jag1"),label = T,max.cutoff = 2)#Macrophage-4
FeaturePlot(scRNA, features = c("C1qa", "Hexb", "Tmem119"),label = T,max.cutoff = 2)#Microglia-0,6
FeaturePlot(scRNA, features = c("H2-Ab1","H2-Aa"),label = T,max.cutoff = 2)#Monocyte DC-4
FeaturePlot(scRNA, features = c("Vtn", "Rgs5", "Acta2"),label = T,max.cutoff = 2)#SMC-1,3,8
FeaturePlot(scRNA, features = c("Tubb2b","Sox11","Meg3"),label = T,max.cutoff = 2)#Neuron-12
FeaturePlot(scRNA, features = c("S100a9"),label = T,max.cutoff = 2)#Neutrophil-11
FeaturePlot(scRNA, features = c("Plp1", "Apod","Cldn11"),label = T,max.cutoff = 2)#Oligodendrocytes-9
FeaturePlot(scRNA, features = c("Ptprz1","Pdgfra","Lhfpl3"),label = T,max.cutoff = 2)#OPC-10
#BM
FeaturePlot(scRNA, features = c("Car2","Tfrc"),label = T,max.cutoff = 2)#MP1--8,13
FeaturePlot(scRNA, features = c("Prtn3","Ctsg"),label = T,max.cutoff = 2)#MP2--10

FeaturePlot(scRNA, features = c("Fcnb","Hmgn2"),label = T,max.cutoff = 2)#N1-4
FeaturePlot(scRNA, features = c("Cebpe1","Ltf"),label = T,max.cutoff = 2)#N2-0,1,3
FeaturePlot(scRNA, features = c("Ltf1","Ngp2"),label = T,max.cutoff = 2)#N3-6
FeaturePlot(scRNA, features = c("Mmp8","Ngp2"),label = T,max.cutoff = 2)#N3-6

FeaturePlot(scRNA, features = c("Lgals11","Crip1"),label = T,max.cutoff = 2)#M1-2,5,11

FeaturePlot(scRNA, features = c("Ccl5","Pdcd4"),label = T,max.cutoff = 2)#LPC-9

FeaturePlot(scRNA, features = c("Cd74","H2-Ab1"),label = T,max.cutoff = 2)#B_cell--7,12


#kidney
Idents(scRNA) <- "integrated_snn_res.0.3"
FeaturePlot(scRNA, features = c("Slc34a1","Lrp2", "Aqp1"),label = T,max.cutoff = 2)#Proximal Tubule Cells--1,9
FeaturePlot(scRNA, features = c("Slc12a1","Umod"),label = T,max.cutoff = 2)#Distal Tubule Cells--2,7,14
FeaturePlot(scRNA, features = c("Aqp2","Hsd11b2"),label = T,max.cutoff = 2)#Collecting Duct Cells-11
FeaturePlot(scRNA, features = c("Pecam1","Cdh5","Emcn"),label = T,max.cutoff = 2)#Endothelial Cells-6,10
FeaturePlot(scRNA, features = c("Adgre1","C1qa","Itgam"),label = T,max.cutoff = 2)#Macrophages-3,5,19,13
FeaturePlot(scRNA, features = c("Cd3e","Cd4","Cd8a"),label = T,max.cutoff = 2)#T Cells-8
FeaturePlot(scRNA, features = c("Cd19","Cd79a","Ms4a1"),label = T,max.cutoff = 2)#B Cells-15
FeaturePlot(scRNA, features = c("Ncr1","Klrb1c","Gzma"),label = T,max.cutoff = 2)#NK cells-17
FeaturePlot(scRNA, features = c("S100a9"),label = T,max.cutoff = 2)#Neutrophils--0,4
FeaturePlot(scRNA, features = c("Itgax"),label = T,max.cutoff = 2)#Dendritic Cells--12
FeaturePlot(scRNA, features = c("Trem1"),label = T,max.cutoff = 2)#Monocytes--13
FeaturePlot(scRNA, features = c("ACta2","Pdgfrb"),label = T,max.cutoff = 2)#Stromal Cells--20
FeaturePlot(scRNA, features = c("Car2","Atp6v1a"),label = T,max.cutoff = 2)#Intercalated Cells--18





scRNA <- subset(scRNA,ident = c(16,21),invert = T)

scRNA@meta.data$celltype <- factor(scRNA@meta.data$integrated_snn_res.0.3,
                                   labels = c("Neutrophils",
                                              "Proximal Tubule Cells","Distal Tubule Cells","Macrophages","Neutrophils",
                                              "Macrophages","Endothelial Cells","Distal Tubule Cells","T Cells",
                                              "Proximal Tubule Cells","Endothelial Cells","Collecting Duct Cells",
                                              "Dendritic Cells","Macrophages","Distal Tubule Cells","B Cells",
                                              "NK cells","Intercalated Cells","Macrophages","Stromal Cells"))
Idents(scRNA) <- "celltype"
DefaultAssay(scRNA) <- "integrated"
DimPlot(scRNA,label = T)
FeaturePlot(scRNA,features = "Trem2")
scRNA_marker <- FindAllMarkers(scRNA,logfc.threshold = 0.5)
write.csv(scRNA_marker,"GSE227075_下肢/marker/marker_celltype.csv")
saveRDS(scRNA,"GSE180420_肾脏2/RDS/kidney_celltype.rds")

DimPlot(scRNA,label = F,group.by = "orig.ident")

#ratio
scRNA <- readRDS("GSE180420_肾脏2/RDS/kidney_mac_celltype.rds")
Idents(scRNA) <- "celltype"
DimPlot(scRNA,label = T)
head(scRNA)
Idents(scRNA) <- "timepoint"
as.data.frame(prop.table(table(scRNA@meta.data[,"celltype"],Idents(scRNA)), margin = 2))-> pdf -> td
library(tidyverse)
library(RColorBrewer)
td$value <- round(td$Freq * 100)
p <- ggplot(td, aes(x=Var2, y=Freq, fill=Var1)) +
  geom_bar(stat="identity", position="stack", aes(fill=Var1)) +
  scale_y_continuous(labels = scales::percent,expand=c(0.01,0.01)) +
  #geom_text(aes(label= Var1), position=position_fill(vjust=0.5))+
  labs(x="sample",y="Cells Ratio")+
  scale_fill_manual(values= brewer.pal(n=12,name = "Paired"))+
  guides(fill = guide_legend(keywidth = 4, keyheight = 1,ncol=1,title = 'celltype',label.position = "bottom"))+
  theme(legend.position = "top")+
  theme_bw()+
  theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p

DimPlot(scRNA,split.by = "timepoint",group.by = "celltype")
#热图
marker <- read.csv("GSE210159_心脏2/marker/marker_mac.csv")
marker <- subset(marker, marker$p_val_adj < 0.05)
top5 <- marker %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)

DefaultAssay(scRNA) <- "integrated"
Idents(scRNA) <- "celltype"
DoHeatmap(scRNA, features = top5$gene)+ scale_fill_gradientn(colors = c("#0099CC","white","#880000"))
FeaturePlot(scRNA,features = "Jag1",label = T)

####mac recluster####
rm(list = ls())
setwd("D:/analysis/公共数据/")
scRNA <- readRDS("GSE227075_下肢/RDS文件/climb_celltype.rds")
Idents(scRNA) <- "celltype"
scRNA$orig.ident <- factor(scRNA@meta.data$orig.ident,
                           levels = c("sham","MI_3D","MI_7D","MI_14D"))
DimPlot(scRNA,label = T)
DimPlot(scRNA,label = F,group.by = "sample")
FeaturePlot(scRNA,features = "Nrg1",label = T)
##提取细胞子集
scRNAsub_EC <- subset(scRNA, idents = "Macrophages")
DefaultAssay(scRNAsub_EC) <- "RNA"
scRNAsub_EC <- NormalizeData(scRNAsub_EC,normalization.method = "LogNormalize",scale.factor = 10000)
scRNAsub_EC <- FindVariableFeatures(scRNAsub_EC, selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(., vars.to.regress = "percent.mt") %>%
  RunPCA() %>%
  RunUMAP(dims = 1:20) %>%
  FindNeighbors(dims = 1:20) %>%
  FindClusters(resolution = c(0.2,0.3,0.4))
Idents(scRNAsub_EC) <- "RNA_snn_res.0.3"
scRNAsub_EC$orig.ident <- factor(scRNAsub_EC@meta.data$orig.ident,
                                 levels = c("sham","MI_3D","MI_7D","MI_14D"))
DimPlot(scRNAsub_EC,label = T,split.by = "sample")
FeaturePlot(scRNAsub_EC,features = "Cd68",label = T)
FeaturePlot(scRNAsub_EC,features = "Trem2",label = T)#0,2
FeaturePlot(scRNAsub_EC,features = "Ltb4r1",label = T)#1
FeaturePlot(scRNAsub_EC,features = "Ifit2",label = T)#9
FeaturePlot(scRNAsub_EC,features = "Lyve1",label = T)#
FeaturePlot(scRNAsub_EC,features = "H2-Eb1",label = T)#3,6
FeaturePlot(scRNAsub_EC,features = "Mki67",label = T)#5
FeaturePlot(scRNAsub_EC,features = "Ly6c2",label = T)#4
FeaturePlot(scRNAsub_EC,features = "S100a9",label = T)#7
FeaturePlot(scRNAsub_EC,features = "Cd209a",label = T)#8

scRNAsub_EC <- subset(scRNAsub_EC,ident = c(8),invert = T)
scRNAsub_EC$celltype <- factor(scRNAsub_EC$RNA_snn_res.0.3,
                               labels = c("Trem2_mac","BLT1_mac","Trem2_mac","MHCII_mac","Ly6c2_mono",
                                          "Ki67_mac","MHCII_mac","S100a9_mono","mono","Ifit2_mac"))
Idents(scRNAsub_EC) <- "celltype"
DimPlot(scRNAsub_EC,label = T)
DimPlot(scRNAsub_EC,label = F,split.by = "sample")
saveRDS(scRNAsub_EC,"GSE180420_肾脏2/RDS/kidney_mac_celltype.rds")

###计算Trem2基因评分
scRNAsub_EC <- readRDS("")
marker <- read.csv("D:/课题数据/ST/analysis/IR五个时间点/marker/mac_marker_celltype_all.csv")
table(marker$cluster)
marker <- subset(marker,marker$cluster == "Trem2_mac" & marker$avg_log2FC > 0 & marker$p_val_adj < 0.05)
gene <- list(marker$gene[1:100])
DefaultAssay(scRNAsub_EC) <- "RNA"
scRNA <- AddModuleScore(scRNAsub_EC,features = gene,name = "Trem2_mac",assay = "RNA")

VlnPlot(scRNA,features = "Trem2_mac1",pt.size = 0.1)

####计算数据相关性####
library(ggpmisc)
library(ggpubr)
library(dplyr)
setwd("D:/analysis/公共数据/")
refdata <- readRDS("D:/课题数据/ST/analysis/IR五个时间点/rds文件/MAC_select_combined.rds")
marker <- read.csv("D:/课题数据/ST/analysis/IR五个时间点/marker/mac marker/mac_marker_celltype_all.csv")
marker <- subset(marker, marker$avg_log2FC > 0 & marker$p_val_adj < 0.05 & marker$cluster == "Trem2_mac")
gene <- marker$gene
gene_exp1 <- data.frame(AverageExpression(refdata,features = gene,group.by = "celltype",assays = "RNA"))
colnames(gene_exp1) <- paste0("Heart_IR_",colnames(gene_exp1))
#gene_exp1$gene <- rownames(gene_exp1)
head(gene_exp1)

scRNA <- readRDS("D:/课题数据/ST/analysis/公共数据/GSE210159_心脏2/RDS/heart_mac_celltype.rds")
gene_exp2 <- data.frame(AverageExpression(scRNA,features = gene,group.by = "celltype",assays = "RNA"))
colnames(gene_exp2) <- paste0("heart_MI",colnames(gene_exp2))
#gene_exp2$gene <- rownames(gene_exp2)
head(gene_exp2)

scRNA <- readRDS("D:/课题数据/ST/analysis/公共数据/GSE180420_肾脏2/RDS/kidney_mac_celltype.rds")
gene_exp3 <- data.frame(AverageExpression(scRNA,features = gene,group.by = "celltype",assays = "RNA"))
colnames(gene_exp3) <- paste0("kidney_",colnames(gene_exp3))
#gene_exp3$gene <- rownames(gene_exp3)
head(gene_exp3)

scRNA <- readRDS("D:/课题数据/ST/analysis/公共数据/GSE227075_下肢/RDS文件/climb_mac_celltype.rds")
gene_exp4 <- data.frame(AverageExpression(scRNA,features = gene,group.by = "celltype",assays = "RNA"))
colnames(gene_exp4) <- paste0("climb_",colnames(gene_exp4))
#gene_exp4$gene <- rownames(gene_exp4)
head(gene_exp4)

gene_exp_all <- Reduce(function(x,y) merge(x,y,by = "gene"),list(gene_exp,gene_exp2,gene_exp3,gene_exp4))
rownames(gene_exp_all) <- gene_exp_all$gene
gene_exp_all$gene <- NULL
#ALL
pheatmap(cor(gene_exp_all),
         scale = "none",
         cluster_rows = F,
         cluster_cols = F,
         gaps_row = c(8,15,23),
         gaps_col = c(8,15,23),
         display_numbers = T)
#Trem2
pheatmap(cor(gene_exp_all[,c(1,9,16,25)]),
         scale = "none",
         cluster_rows = F,
         cluster_cols = F,
         #gaps_row = c(8,15,23),
         #gaps_col = c(8,15,23),
         display_numbers = T)
#BLT1
pheatmap(cor(gene_exp_all[,c(5,10,17,26)]),
         scale = "none",
         cluster_rows = F,
         cluster_cols = F,
         #gaps_row = c(8,15,23),
         #gaps_col = c(8,15,23),
         display_numbers = T)


gene_exp1 <- subset(gene_exp1,rownames(gene_exp1) %in% rownames(gene_exp4))
data_plot <- cbind(gene_exp1,gene_exp4)
cor=round(cor(data_plot$Heart_IR_RNA.Trem2_mac,
              data_plot$climb_RNA.Trem2_mac),
          2)
ggplot(log1p(data_plot), aes(x = Heart_IR_RNA.Trem2_mac, y = climb_RNA.Trem2_mac)) +
  geom_point(color="#6baed6")+
  geom_smooth(method = "lm", formula = y~x, color = "#756bb1", fill = "#cbc9e2")+ 
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  stat_cor(size = 8)

####基因表达####
setwd("H:/Usb Drivе/课题数据/ST/analysis/公共数据/")
scRNA <- readRDS("GSE210159_心脏2/RDS/heart_mac_celltype.rds")
head(scRNA)
VlnPlot(scRNA,"Pf4")
colnames(scRNA@meta.data)[14] <- "group"
table(scRNA@meta.data$group)
scRNA@meta.data$group <- factor(scRNA@meta.data$group,
                                labels = c("Heart_sham","Heart_1D","Heart_3D"))
scRNA$celltype <- factor(scRNA$celltype,levels = c("Trem2_mac","MHCII_mac","Lyve1_mac","BLT1_mac","Ifit2_mac","Ki67_mac"))
Idents(scRNA) <- "celltype"
DimPlot(scRNA,label = F,pt.size = 1)
saveRDS(scRNA,"GSE210159_心脏2/RDS/heart_mac_celltype.rds")
scRNA <- subset(scRNA,ident = "mono",invert = T)
DimPlot(scRNA,label = T)
DefaultAssay(scRNA) <- "RNA"
genelist <- c("Trem2","H2-Eb1","Lyve1","Ltb4r1","Ifit2","Mki67")
VlnPlot(scRNA,features = genelist,pt.size = 0,stack = T) +
  stat_summary(fun.y=median, geom="point",size=3, color="black")
p1 <- FeaturePlot(scRNA,features = "Jag1",label = T)
p2 <- VlnPlot(scRNA,features = "Jag1",pt.size = 0.1)
p3 <- DimPlot(scRNA,label = T)
p1 + p2 + p3

####比例变化####
setwd("D:/课题数据/ST/analysis/公共数据/")
scRNA <- readRDS("GSE227075_下肢/RDS文件/climb_mac_celltype.rds")
DimPlot(scRNA,label = T)
#单个细胞类型
# 1. 提取元数据并筛选目标细胞类型
metadata <- scRNA@meta.data
target_type <- "Trem2_mac"  # 修改为目标细胞类型名称

# 2. 统计每个分组中的细胞数量
cell_counts <- metadata %>%
  group_by(group) %>%        # 按分组统计
  summarise(
    total_cells = n(),      # 总细胞数
    target_cells = sum(celltype == target_type)  # 目标细胞数
  ) %>%
  mutate(
    proportion = target_cells / total_cells * 100  # 计算百分比
  )

# 3. 确保分组顺序正确（如果是时间序列或有序分类）
cell_counts$group <- factor(cell_counts$group, levels = cell_counts$group)  # 示例分组顺序

# 4. 绘制折线图
ggplot(cell_counts, aes(x = group, y = proportion, group = 1)) +
  geom_line(color = "#2c7bb6", linewidth = 1.5) +       # 蓝色折线
  geom_point(color = "#d7191c", size = 4) +            # 红色点
  labs(
    x = "Experimental Group",
    y = "Cell Proportion (%)",
    title = paste("Proportion of", target_type, "Across Groups")
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.title = element_text(size = 12)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))  # 调整Y轴范围
#全部细胞类型
scRNA <- readRDS("GSE210159_心脏2/RDS/heart_mac_celltype.rds")
DimPlot(scRNA,label = T)
# 1. 提取元数据并统计各细胞类型的比例
cell_prop <- scRNA@meta.data %>%
  group_by(group, celltype) %>%        # 按分组和细胞类型统计
  summarise(count = n(), .groups = "drop") %>%  # 计算每个组内各类型的细胞数
  group_by(group) %>%                  # 按组计算比例
  mutate(proportion = count / sum(count) * 100) %>%
  ungroup()

# 2. 确保分组顺序正确（如果是时间序列或有序分类）
cell_prop$group <- factor(
  cell_prop$group, 
  levels = unique(cell_prop$group)  # 按实际分组顺序调整
)
# 3. 自定义颜色方案（按细胞类型数量自动生成）
celltypes <- unique(cell_prop$celltype)
color_palette <- scales::hue_pal()(length(celltypes))  # 使用 ggplot 默认色盘

# 4. 绘制所有细胞类型的折线图
ggplot(cell_prop, aes(x = group, y = proportion, color = celltype, group = celltype)) +
  geom_line(linewidth = 1.2) +                      # 折线
  geom_point(size = 3) +                            # 数据点
  labs(
    x = "Experimental Group",
    y = "Cell Proportion (%)",
    color = "Cell Type",
    title = "Proportion of All Cell Types Across Groups"
  ) +
  scale_color_manual(values = color_palette) +       # 自定义颜色
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    legend.key.height = unit(0.5, "cm"),            # 紧凑图例
    legend.text = element_text(size = 8)
  ) 
##Ro/e使用折线图进行可视化
scRNA <- readRDS("GSE227075_下肢/RDS文件/climb_mac_celltype.rds")
metadata=scRNA@meta.data
head(scRNA)
summary <- table(metadata[,c('celltype','group')])

roe <- as.data.frame(t(ROIE(summary)))
roe$group <- rownames(roe)
roe$group <- factor(roe$group,levels = unique(roe$group))
ggplot(roe, aes(x = roe$group, y = Trem2_mac, group = 1)) +
  geom_line(color = "#2c7bb6", linewidth = 1.5) +       # 蓝色折线
  geom_point(color = "#d7191c", size = 4) +            # 红色点
  labs(
    x = "Group",
    y = "Ro/e",
    title = paste("Ro/e of Trem2_mac Across Groups")
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.title = element_text(size = 12)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

####cellchat####
#数据处理
scRNA <- readRDS("GSE180420_肾脏2/RDS/kidney_celltype.rds")
DimPlot(scRNA,label = T,group.by = "celltype")

scRNAsub <- readRDS("GSE180420_肾脏2/RDS/kidney_mac_celltype.rds")
DimPlot(scRNAsub,label = T,group.by = "celltype")

df1 <- data.frame(scRNA@meta.data) %>% rownames_to_column()
head(df1)
df1 <- df1[,c(1,11)]
colnames(df1) <- c("key_column","value_column")
df2 <- data.frame(scRNAsub@meta.data) %>% rownames_to_column()
df2 <- df2[,c(1,11)]
colnames(df2) <- c("key_column","value_column")
head(df2)

# 使用 left_join 结合 mutate 和 case_when 来替换值
result <- df1 %>%
  left_join(df2, by = "key_column", suffix = c("_df1", "_df2")) %>%
  mutate(
    value_column = case_when(
      !is.na(value_column_df2) ~ value_column_df2,  # 如果 df2 中存在对应值，则替换
      TRUE ~ value_column_df1                        # 否则保留原 df1 的值
    )
  ) %>%
  select(key_column, value_column)

# 输出结果
print(result)
scRNA$cellchat <- result$value_column
DimPlot(scRNA,group.by = "cellchat",label = T)

#cellchat
library(CellChat)
##创建cellchat对象
Idents(scRNA) <- "cellchat"
scRNA@meta.data$cellchat <- factor(scRNA@meta.data$cellchat,levels = unique(scRNA@meta.data$cellchat))
data.input  <- scRNA@assays$RNA@data
identity = data.frame(group =scRNA$cellchat, row.names = names(scRNA$cellchat)) # create a dataframe consisting of the cell labels
unique(identity$group) # check the cell labels
cellchat <- createCellChat(data.input)
#把metadata信息加到CellChat对象
cellchat <- addMeta(cellchat, meta = identity, meta.name = "labels")
cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group
##导入配受体数据库
CellChatDB <- CellChatDB.mouse
#从特定的方向来刻画细胞间相互作用
#可选择得方向
unique(CellChatDB$interaction$annotation)
#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling for cell-cell communication analysis
cellchat@DB <- CellChatDB # set the used database in the object

##预处理,提取过表达得配体-受体对和相互作用
cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
#识别过表达的配体或受体，然后将基因表达数据投射到蛋白-蛋白相互作用(PPI)网络上
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.mouse)  

##相互作用推断
#为每个相互作用分配一个概率值并进行置换检验来推断生物意义上的细胞-细胞通信
#每个配体-受体对的细胞间通信网络和每个信号通路分别存储在“net”和“netP”槽中
cellchat <- computeCommunProb(cellchat,type = "truncatedMean",trim = 0.1)
#在信号通路级别推断细胞-细胞通信
cellchat <- computeCommunProbPathway(cellchat)
#计算整合的细胞通信网络
cellchat <- aggregateNet(cellchat)
#数据保存
df.net <- subsetCommunication(cellchat)
write.csv(df.net,"GSE180420_肾脏2/cellchat分析/L_R.csv")
saveRDS(cellchat,"GSE180420_肾脏2/RDS/climb_cellchat.rds")
unique(scRNA@meta.data$cellchat)
netVisual_bubble(cellchat, sources.use = c(14), targets.use = c(6,11), remove.isolate = FALSE)

#cellchat统计
library(RColorBrewer)
data <- read.csv("GSE227075_下肢/cellchat/L_R.csv")
head(data)
table(data$target)
celltype_climb <- c("Cyc_cell","EC","Epi_cell","FAPs","MuSCs")
datasub <- subset(data,data$source == "Trem2_mac" & data$target %in% celltype)
sum <- data.frame(table(datasub$target))
head(sum)
sum$Var1 <- factor(sum$Var1,levels = celltype)
ggplot(sum, mapping=aes(x = Var1, y = Freq,fill = Var1))+
  geom_bar(stat="identity")+
  labs(x = "celltype", y = "receptor number",title = "climb Trem2_mac-other") + 
  theme_bw()+
  scale_fill_manual(values = brewer.pal(11,"Paired")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,face = "bold",size = 12))
####功能富集分析(GO)####
rm(list = ls())
library(org.Mm.eg.db)
library(clusterProfiler)
library(dplyr)
#marker计算
setwd("D:/analysis/公共数据/")
scRNA <- readRDS("GSE227075_下肢/RDS文件/climb_mac_celltype.rds")
Idents(scRNA) <- "timepoint"
scRNA <- subset(scRNA,ident = "Day3")
Idents(scRNA) <- "celltype"
DimPlot(scRNA,label = T)
DefaultAssay(scRNA) <- "integrated"
marker <- FindAllMarkers(scRNA)
marker <- subset(marker,marker$p_val_adj < 0.05)
write.csv(marker,"GSE227075_下肢/marker/marker_mac_72h.csv")
#marker基因比较
library(VennDiagram)
marker_heart_IR <- read.csv("心脏IR/marker/mac_marker_celltype.csv")
marker_heart_IR_Trem2 <- subset(marker_heart_IR,marker_heart_IR$cluster == "Trem2_mac")$gene
marker_climb <- read.csv("GSE227075_下肢/marker/marker_mac.csv")
marker_climb_Trem2 <- subset(marker_climb,marke_climb$cluster == "Trem2_mac")$gene
marker_heart_MI <- read.csv("GSE210159_心脏2/marker/marker_mac.csv")
marker_heart_MI_Trem2 <- subset(marker_heart_MI,marker_heart_MI$cluster == "Trem2_mac")$gene

intersection_trem2 <- Reduce(intersect,list(marker_heart_IR_Trem2[1:100],marker_climb_Trem2[1:100]))

draw.pairwise.venn(
  area1 = length(marker_heart_IR_Trem2),
  area2 = length(marker_heart_MI_Trem2),
  cross.area = length(intersection_trem2),
  category = c("heart_IR_Trem2", "heart_MI_Trem2"),
  fill = c("red", "blue"),
  lty = "blank",
  cex = 1,
  cat.cex = 1,
  cat.col = c("red", "blue")
)
dev.off()
#GO分析
empty_df <- data.frame()
cluster_name <- data.frame(table(marker$cluster))$Var1
for(i in cluster_name){
  genelist <- subset(marker,marker$cluster == i & marker$p_val_adj < 0.05 & marker$avg_log2FC > 0 )
  eg <- bitr(genelist$gene, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Mm.eg.db")
  id=na.omit(eg$ENTREZID)
  BP <- enrichGO(id, "org.Mm.eg.db", keyType = "ENTREZID",ont = 'BP',pvalueCutoff  = 0.05,pAdjustMethod = "BH",  qvalueCutoff  = 0.1, readable=T)
  BP@result$celltype <- rep(paste0("climb_",i),length(rownames(BP@result)))
  BP@result <- subset(BP@result, BP@result$p.adjust < 0.05)
  empty_df <- rbind(empty_df,BP@result)
  write.csv(BP@result, file=paste0("GSE227075_下肢/GO分析/mac_72h_",i,"_GO_BP.csv"))
}
write.csv(empty_df,"GSE227075_下肢/GO分析/mac_72h_all_mac_GO_BP.csv")

####功能富集分析(GSEA-marker)####
library(clusterProfiler)
library(org.Mm.eg.db)
library(stringr)
library(ggplot2)
library(enrichplot)
library(msigdbr)
library(tidyr)
library(tibble)
rm(list = ls())
marker <- read.csv("GSE227075_下肢/marker/marker_mac_72h.csv")
cluster_name <- data.frame(table(marker$cluster))$Var1
all_gsea <- data.frame()
for(i in 1:6){
  marker_Trem2 <- subset(marker,marker$cluster == cluster_name[i] & marker$p_val_adj < 0.05)
  marker_Trem2 <- marker_Trem2[order(-marker_Trem2$avg_log2FC),]
  marker_Trem2 <- dplyr::rename(marker_Trem2,SYMBOL = gene)
  gene <- str_trim(marker_Trem2$SYMBOL,"both")
  #ID转换
  gene=bitr(gene,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Mm.eg.db") 
  ## 去重
  gene <- dplyr::distinct(gene,SYMBOL,.keep_all=TRUE)
  gene_df <- merge(marker_Trem2,gene,by="SYMBOL")
  #定义基因列表，对logFC进行从高到低排序
  geneList <- gene_df$avg_log2FC #第二列可以是folodchange，也可以是logFC
  names(geneList)=gene_df$ENTREZID #使用转换好的ID
  geneList=sort(geneList,decreasing = T) #从高到低排序
  gse.GO <- gseGO(
    geneList, #geneList
    ont = "BP",  # 可选"BP"、"MF"和"CC"或"ALL"
    OrgDb = org.Mm.eg.db, #人 注释基因
    keyType = "ENTREZID",
    pvalueCutoff = 0.5,
    pAdjustMethod = "BH",#p值校正方法
  )
  write.csv(gse.GO@result,paste0("GSE227075_下肢/GSEA分析72H/",cluster_name[i],"_GSEA.csv"))
  gse.GO@result$cluster <- paste(cluster_name[i])
  all_gsea <- rbind(all_gsea,gse.GO@result)
}
write.csv(all_gsea,"GSE227075_下肢/GSEA分析72H/all_heart_mac_GSEA.csv")
#可视化
go_data <- read.csv("GSE227075_下肢/GSEA分析72H/Trem2_mac_GSEA.csv")
term <- read.csv("GSE227075_下肢/GSEA分析72H/term.csv")
filtered_data <- subset(go_data,go_data$Description %in% term$term)
filtered_data$Description <- factor(filtered_data$Description,levels = rev(filtered_data$Description))
ggplot(filtered_data, aes(
  x = reorder(Description, enrichmentScore),  # 按 enrichmentScore 排序通路名称
  y = enrichmentScore)) +
  geom_bar(stat = "identity",fill = "#c85e62") +
  labs(
    x = "Pathway",
    y = "enrichmentScore",
    title = "Significantly Enriched Pathways of Trem2_mac (NES > 1)"
  ) +
  coord_flip() +  # 翻转坐标轴使长名称可读
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    axis.text.y = element_text(size = 16,face = "bold"),
    legend.position = "top"
  )



#####功能富集分析(GSEA-expres level)####
library(VennDiagram)
IR_mac <- readRDS("D:/课题数据/ST/analysis/IR五个时间点/rds文件/MAC_select_combined.rds")
expression1 <- data.frame(AverageExpression(IR_mac,group.by = "celltype")$RNA)
expression1 <- expression[order(-expression1$Trem2_mac),]
climb_mac <- readRDS("GSE227075_下肢/RDS文件/climb_mac_celltype.rds")
expression2 <- data.frame(AverageExpression(climb_mac,group.by = "celltype")$RNA)
expression2 <- expression[order(-expression2$Trem2_mac),]

intersection_trem2 <- Reduce(intersect,list(rownames(expression1)[1:1000],rownames(expression2)[1:1000]))
draw.pairwise.venn(
  area1 = length(rownames(expression1)[1:1000]),
  area2 = length(rownames(expression2)[1:1000]),
  cross.area = length(intersection_trem2),
  category = c("heart_IR_Trem2", "heart_MI_Trem2"),
  fill = c("red", "blue"),
  lty = "blank",
  cex = 1,
  cat.cex = 1,
  cat.col = c("red", "blue")
)
dev.off()
#计算差异基因交集
library(tidyr)
library(dplyr)
heart <- read.csv("心脏IR/GSEA分析/mac_72h_Trem2_mac_GSEA.csv")
heart <- heart[,c(3,5,13)]
heart2 <- read.csv("GSE210159_心脏2/GSEA分析/mac_Trem2_mac_GSEA.csv")
heart2 <- heart2[,c(3,5,13)]
brain <- read.csv("GSE210986_脑/GSEA分析/mac_Trem2_mac_GSEA.csv")
brain <- brain[,c(3,5,13)]
climb <- read.csv("GSE227075_下肢/GSEA分析/mac_72h_Trem2_mac_GSEA.csv")
climb <- climb[,c(3,5,13)]
kidney <- read.csv("GSE267242_肾脏/GSEA分析/mac_Trem2_mac_GSEA.csv")
kidney <- kidney[,c(3,5,13)]
GSEAdata <- rbind(heart,climb)
GSEAdata_plot <- pivot_wider(GSEAdata,names_from = celltype,values_from = enrichmentScore,values_fill = 0)
GSEAdata_plot <- GSEAdata_plot %>% column_to_rownames(var = "Description")
pheatmap(GSEAdata_plot,
         show_rownames = F)


heart <- read.csv("心脏IR/GO分析/mac_Trem2_mac_GO_BP.csv")
head(heart)

heart$p.adjust <- -log10(heart$p.adjust)
heart <- heart[1:20,]
heart$Description <- factor(heart$Description,levels = rev(heart$Description))
# 绘制气泡图
ggplot(heart, aes(x = p.adjust, y = Description, size = Count, fill = p.adjust)) +
  geom_point(shape = 21, color = "black", stroke = 0.5) +
  scale_size(range = c(3, 6)) +  # 调整气泡的大小范围
  scale_fill_gradient(low = "lightblue", high = "red") +  # 设置颜色渐变
  theme_minimal() +
  labs(title = "heart Trem2 Mac GO analysis",
       x = "-log10(p-value)",
       y = "GO Terms",
       size = "Gene Count",
       fill = "-log10(p-value)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(face = "bold"))


####细胞群相似性计算#####
library(Seurat)
library(pheatmap)
setwd("D:/analysis/公共数据/")
seurat1 <- readRDS("D:/课题数据/ST/analysis/IR五个时间点/rds文件/MAC_select_combined.rds")
DimPlot(seurat1)
seurat2 <- readRDS("GSE180420_肾脏2/RDS/kidney_mac_celltype2.rds")
DimPlot(seurat2)
# 假设 data1 和 data2 是两个 Seurat 对象
# Step 1: 确保两个数据集有共同的基因
common_genes <- intersect(rownames(seurat1), rownames(seurat2))

# 只保留共同基因
seurat1 <- subset(seurat1, features = common_genes)
seurat2 <- subset(seurat2, features = common_genes)

# Step 2: 按群体计算平均表达值
# Step 1: 确认两个对象中存在的共同注释类别
# 假设群体注释存储在 meta.data 中的 "cluster" 列
common_celltypes <- intersect(unique(seurat1$celltype), unique(seurat2$celltype))
seurat1 <- subset(seurat1, cells = colnames(seurat1)[seurat1$celltype %in% common_celltypes])
seurat2 <- subset(seurat2, cells = colnames(seurat2)[seurat2$celltype %in% common_celltypes])

# 计算每个群体的平均表达值
average_expression1 <- AverageExpression(seurat1, group.by = "celltype", assay = "RNA")$RNA
average_expression2 <- AverageExpression(seurat2, group.by = "celltype", assay = "RNA")$RNA

# Step 3: 计算相关性矩阵
correlation_matrix <- cor(average_expression1, average_expression2, method = "pearson")

# Step 4: 可视化相关性矩阵
pheatmap(correlation_matrix, 
         cluster_rows = F, 
         cluster_cols = F, 
         display_numbers = TRUE, 
         color = colorRampPalette(c("blue", "white", "red"))(50), 
         main = "Correlation between Cell Groups",
         fontsize = 15)
#降维
seurat1 <- readRDS("D:/课题数据/ST/analysis/IR五个时间点/rds文件/MAC_select_combined.rds")
DimPlot(seurat1)
seurat1$dataorig <- c(rep("H_IR",length(colnames(seurat1))))
seurat2 <- readRDS("GSE210159_心脏2/RDS/heart_mac_celltype.rds")
DimPlot(seurat2)
seurat2$dataorig <- c(rep("H_MI",length(colnames(seurat2))))
seurat3 <- readRDS("GSE180420_肾脏2/RDS/kidney_mac_celltype2.rds")
DimPlot(seurat3)
seurat3$dataorig <- c(rep("K_IR",length(colnames(seurat3))))
seurat4 <- readRDS("GSE227075_下肢/RDS文件/climb_mac_celltype.rds")
DimPlot(seurat4)
seurat4$dataorig <- c(rep("C_IR",length(colnames(seurat4))))
# 假设 data1 和 data2 是两个 Seurat 对象
# Step 1: 确保两个数据集有共同的基因
list <- list(rownames(seurat1), rownames(seurat2), rownames(seurat3), rownames(seurat4))
common_genes <- reduce(list,intersect)

# 只保留共同基因
seurat1 <- subset(seurat1, features = common_genes)
DefaultAssay(seurat1) <- "RNA"
seurat2 <- subset(seurat2, features = common_genes)
DefaultAssay(seurat2) <- "RNA"
seurat3 <- subset(seurat3, features = common_genes)
DefaultAssay(seurat3) <- "RNA"
seurat4 <- subset(seurat4, features = common_genes)
DefaultAssay(seurat4) <- "RNA"

split_seurat <- list(seurat1,seurat2,seurat3,seurat4)
library(future)
library(harmony)
options(future.globals.maxSize = 8000 * 1024^2)
# normalize and identify variable features for each dataset independently
ifnb.list <- lapply(X = split_seurat, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# Step 4: 识别整合锚点
# 使用 Seurat 的整合方法
anchors <- FindIntegrationAnchors(object.list = ifnb.list, 
                                  dims = 1:30)  # dims 表示用于整合的 PCA 维度

# Step 5: 整合数据
# 整合后的对象会消除批次效应
integrated_seurat <- IntegrateData(anchorset = anchors, dims = 1:30)

# Step 6: 后续分析（标准化、降维等）
# 设置默认 Assay 为整合后的数据
DefaultAssay(integrated_seurat) <- "integrated"

# 数据缩放
integrated_seurat <- ScaleData(integrated_seurat)

# PCA 降维
integrated_seurat <- RunPCA(integrated_seurat)

# UMAP 降维
integrated_seurat <- RunUMAP(integrated_seurat, dims = 1:30)

# Step 5: 可视化降维结果
# 使用保留的原始注释进行可视化
head(integrated_seurat)
integrated_seurat$celltype_new <- paste0(integrated_seurat$dataorig,integrated_seurat$celltype)
DimPlot(integrated_seurat, reduction = "umap", group.by = "celltype", label = T) + scale_color_manual(values = sci_palette)
DimPlot(integrated_seurat, reduction = "umap", group.by = "dataorig", label = F)
saveRDS(integrated_seurat,"integrated.rds")
####figure2####
sci_palette <- c(
  "#D1352B","#D2EBC8","#7DBFA7","#EE934E","#9B5B33","#B383B9","#FCED82","#3C77AF","#AECDE1","#F5CFE4","#8FA4AE","#F5D2A8","#BBDD78"
  )
setwd("D:/课题数据/ST/analysis/公共数据/")
scRNA_all <- readRDS("GSE210159_心脏2/RDS/")
scRNA_all$celltype <- factor(scRNA_all$celltype,levels = c("Trem2_mac","MHCII_mac",
                                                           "Lyve1_mac","Ly6c2_mono","BLT1_mac"
                                                           ,"Ifit2_mac","S100a9_mono","Ki67_mac"))
Idents(scRNA_all) <- "celltype"
DimPlot(scRNA_all,label = T) + scale_color_manual(values = sci_palette)
DimPlot(scRNA_all,label = F,group.by = "dataorig") + scale_color_manual(values = sci_palette)

scRNA <- readRDS("GSE227075_下肢/RDS文件/climb_celltype.rds")
DimPlot(scRNA,label = T) + scale_color_manual(values = sci_palette)
scRNAsub <- readRDS("GSE227075_下肢/RDS文件/climb_mac_celltype.rds")
DimPlot(scRNAsub,label = T) + scale_color_manual(values = sci_palette)
#R o/e
metadata=scRNAsub@meta.data
head(scRNAsub)
summary <- table(metadata[,c('celltype','group')])

roe <- as.data.frame(ROIE(summary))

library(pheatmap)
library(RColorBrewer)
#data <- roe[-6,]
pheatmap(roe, 
         scale = "row",
         display_numbers = TRUE,
         number_color = "black",
         cluster_row = FALSE,
         cluster_col = FALSE,main="Ro/e of C_IR",
         colorRampPalette(c("#8FA4AE","#AECDE1","#D1352B"))(50))

####各亚群GSEA分析####
library(tibble)
allgsea <- read.csv("公共数据/GSE210159_心脏2/GSEA分析/all_heart_mac_GSEA.csv")
allgsea <-  allgsea %>% arrange(desc(allgsea$enrichmentScore))
term <- read.csv("公共数据/GSE210159_心脏2/GSEA分析/term.csv")
term_name <- term$GSEA[1:7]
cluster_name <- term$macsub
gseaplot <- subset(allgsea,allgsea$Description %in% term_name)
gseaplot <- gseaplot[,c(3,5,13)]

GSEAdata_plot <- pivot_wider(gseaplot ,names_from = cluster,values_from = enrichmentScore,values_fill = 0)
GSEAdata_plot <- GSEAdata_plot %>% column_to_rownames(var = "Description")
GSEAdata_plot2 <- GSEAdata_plot[term_name,cluster_name]
pheatmap(GSEAdata_plot2,
         cluster_rows = F,
         cluster_cols = F,
         scale = "none",
         fontsize_row = 16,
         fontsize_col = 16,
         main = "GSEA of macsub in heart_MI")

d# 加载必要的包
library(ggplot2)
library(dplyr)

# 1. 读取数据
allgsea <- read.csv("公共数据/GSE227075_下肢/GSEA分析/all_heart_mac_GSEA.csv")
term <- read.csv("公共数据/GSE227075_下肢/GSEA分析/term.csv")

# 2. 提取需要展示的通路和细胞亚群顺序
term_name <- term$GSEA[1:6]
cluster_name <- unique(term$macsub) # 使用 unique 防止重复项导致 factor 报错

# 3. 过滤并提取长表数据
# 强烈建议：放弃原代码中 [,c(3,5,13)] 这种数字索引，改用列名提取，避免错乱
gseaplot <- allgsea %>%
  filter(Description %in% term_name) %>%
  select(cluster, Description, NES, setSize, p.adjust) # 假设你的列名包含这些，若不同请微调

# 4. 严格固定横轴和纵轴的顺序
# 纵轴 (Y轴): 通路名称，使用 rev() 反转顺序，确保 term_name 的第一个元素出现在图片的顶端
gseaplot$Description <- factor(gseaplot$Description, levels = term_name)

# 横轴 (X轴): 细胞亚群，按照你从 term$macsub 中提取的顺序排列
gseaplot$cluster <- factor(gseaplot$cluster, levels = cluster_name)

# 5. 绘制气泡图
ggplot(gseaplot, aes(x = cluster, y = Description)) +
  # 颜色映射为 NES，大小映射为 setSize
  geom_point(aes(color = NES, size = setSize)) +
  # 颜色梯度：红正蓝负
  scale_color_gradient2(low = "#3182BD", mid = "#F0F0F0", high = "#DE2D26", midpoint = 0) +
  # 调整气泡显示范围
  scale_size_continuous(range = c(3, 9)) +
  # 图表主题美化
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title = element_blank(),
    panel.grid.major = element_line(color = "grey90", linetype = "dashed"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  labs(
    title = "GSEA of macsub in heart_MI",
    color = "NES",
    size = "Gene Set Size"
  )

# 如果需要保存图片，取消下方代码的注释
ggsave("GSEA_climb_macsub_bubble.pdf", width = 8, height = 4,device = cairo_pdf)
