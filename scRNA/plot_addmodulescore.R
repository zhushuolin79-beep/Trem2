library(Seurat)

# trem_2_top_20_de <- c('Spp1', 'Fabp5', 'Lgals3', 'Pf4', 'Ccl7', 'Ctsd', 'Ctsb', 'Gpnmb', 'Lgmn', 'Ctsl', 'Hmox1', 'Trem2', 'Fn1', 'Cstb', 'Cd68', 'Arg1', 'Rgs1', 'Ccl8', 'Plin2', 'Ftl1')
# lyve1_top_20_de <- c('Folr2', 'Lyve1', 'C1qa', 'Cbr2', 'Apoe', 'C1qc', 'C1qb', 'Sepp1', 'F13a1', 'Pltp', 'Mrc1', 'Cd163', 'Lyz2', 'Pf4', 'Ltc4s', 'Trf', 'Maf', 'Fcgrt', 'Mgl2', 'Fcrls')
# mhc_top_20_de <- c('H2-Eb1', 'H2-Aa', 'H2-Ab1', 'Cd74', 'C1qa', 'C1qb', 'C1qc', 'Apoe', 'Mgl2', 'Tmem176b', 'Tmem176a', 'Cd72', 'Trf', 'Cd83', 'H2-DMb1', 'Ccl12', 'Csf1r', 'Aif1', 'Slamf9', 'Cxcl16')

mac_de <- read.csv('mac_top20_de.txt', sep = '\t', header = F)

temp = split(mac_de, mac_de$V1)

ir72 <-readRDS('IR72h.rds')
DefaultAssay(ir72) <- 'Spatial'
Idents(ir72) <- factor(ir72@meta.data$SpaGCN_refined_pred, levels = seq(0,6))

for (item in temp){
  ct = unique(item$V1)
  genes = item$V2
  ir72 <- AddModuleScore(ir72, features = list(genes), name = paste0(ct, '_score_'))
}

ec_de <- read.csv('ec_top20_de.txt', sep = '\t', header = F)

temp = split(ec_de, ec_de$V1)

ir72 <-readRDS('../analysis_7_4/rds_update/IR72h_with_major_frac.rds')
DefaultAssay(ir72) <- 'Spatial'
Idents(ir72) <- factor(ir72@meta.data$SpaGCN_refined_pred, levels = seq(0,6))

for (item in temp){
  ct = unique(item$V1)
  print (ct)
  genes = item$V2
  ir72 <- AddModuleScore(ir72, features = list(genes), name = paste0(ct, '_score_'))
}


ec_for_vln <- read.csv('ec_for_vln.txt', sep = '\t')
ggplot(cf_for_vln, aes(x=Type, y=Score, fill=Type)) + geom_violin(scale = 'width')+ geom_boxplot(width=0.1, outlier.alpha = 0) + facet_wrap(~Cluster) + theme(axis.text.x = element_text(angle=45, hjust = 1))

cf_de <- read.csv('cf_top20_de.txt', sep = '\t', header = F)

temp = split(cf_de, cf_de$V1)

ir72 <-readRDS('../analysis_7_4/rds_update/IR72h_with_major_frac.rds')
DefaultAssay(ir72) <- 'Spatial'
Idents(ir72) <- factor(ir72@meta.data$SpaGCN_refined_pred, levels = seq(0,6))

for (item in temp){
  ct = unique(item$V1)
  print (ct)
  genes = item$V2
  ir72 <- AddModuleScore(ir72, features = list(genes), name = paste0(ct, '_score_'))
}


ec_for_vln <- read.csv('ec_for_vln.txt', sep = '\t')
ggplot(ec_for_vln, aes(x=Type, y=Score, fill=Type)) + 
  geom_violin(scale = 'width')+ 
  geom_boxplot(width=0.1, outlier.alpha = 0) + 
  facet_wrap(~Cluster) + 
  theme(axis.text.x = element_text(angle=45, hjust = 1))
