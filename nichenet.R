rm(list = ls())
library(nichenetr)
library(Seurat)
library(tidyverse)
library(cowplot)

setwd("E:/ST/")

rm(list = ls())
ir72 <- readRDS('analysis/IR五个时间点/rds文件/IR72h_H.rds')
ir72 <- SetIdent(ir72, value = ir72@meta.data$cellchat_Mac)
#Ligand-target model
ligand_target_matrix <- readRDS('Company results/nichenet/ligand_target_matrix.rds')
ligand_target_matrix[1:5,1:5]
#NicheNet’s ligand-receptor data sources
lr_network <- readRDS('Company results/nichenet/lr_network.rds')
weighted_networks <- readRDS('Company results/nichenet/weighted_networks.rds')

weighted_networks_lr = weighted_networks$lr_sig %>% 
  inner_join(lr_network %>% distinct(from,to), by = c("from","to"))
#人转鼠
lr_network = lr_network %>% 
  mutate(from = convert_human_to_mouse_symbols(from), to = convert_human_to_mouse_symbols(to)) %>% 
  drop_na()

colnames(ligand_target_matrix) = ligand_target_matrix %>% 
  colnames() %>% 
  convert_human_to_mouse_symbols()

rownames(ligand_target_matrix) = ligand_target_matrix %>% 
  rownames() %>% 
  convert_human_to_mouse_symbols()
#删除NA值
ligand_target_matrix = ligand_target_matrix %>% 
  .[!is.na(rownames(ligand_target_matrix)), !is.na(colnames(ligand_target_matrix))]

weighted_networks_lr = weighted_networks_lr %>% 
  mutate(from = convert_human_to_mouse_symbols(from), 
         to = convert_human_to_mouse_symbols(to)) %>% drop_na()
ligand_target_matrix[1:5,1:5]

#受体细胞表达基因集
celltype <- c("Bcell","DC","LPC","Neu","Tcell","CF", "CM","EC","SMC","Pericyte")
receiver = celltype
expressed_genes_receiver = get_expressed_genes(receiver, ir72, pct = 0.10, assay_oi = 'RNA')
background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
#mac <- c("BLT1_mac","Ifit2_mac","Trem2_mac","MHCII_mac",
#         "Lyve1_mac","Ki67_mac","S100a9_mono","Ly6c2_mono")
#配体细胞表达基因集
sender_celltypes = mac
list_expressed_genes_sender = sender_celltypes %>% unique() %>% lapply(get_expressed_genes, ir72, 0.10, 'RNA') # lapply to get the expressed genes of every sender cell type separately here
expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()
#受体细胞marker基因
DE_table_receiver = FindMarkers(object = ir72, ident.1 = receiver, min.pct = 0.1) %>% 
  rownames_to_column('gene')
geneset_oi = DE_table_receiver %>% 
  filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25) %>% 
  pull(gene)
#筛选存在于受配体数据库中的marker基因
geneset_oi = geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]

#受配体细胞基因与数据库受配体基因取交集
ligands = lr_network %>% pull(from) %>% unique()
receptors = lr_network %>% pull(to) %>% unique()
expressed_ligands = intersect(ligands,expressed_genes_sender)
expressed_receptors = intersect(receptors,expressed_genes_receiver)

#筛选存在受体的潜在配体
potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()
#预测配体在调节感兴趣基因组表达中的活性
ligand_activities = predict_ligand_activities(geneset = geneset_oi, #受体细胞marker基因
                                              background_expressed_genes = background_expressed_genes, #目标受体基因
                                              ligand_target_matrix = ligand_target_matrix, #参考数据集
                                              potential_ligands = potential_ligands) #潜在的配体基因
ligand_activities = ligand_activities %>% arrange(-pearson) %>% mutate(rank = rank(desc(pearson)))
ligand_activities
#选取top20
best_upstream_ligands = ligand_activities %>% top_n(20, pearson) %>% arrange(-pearson) %>% pull(test_ligand) %>% unique()
DotPlot(ir72, features = best_upstream_ligands %>% rev(), cols = "RdYlBu") + RotatedAxis()

#Active target gene inference
#get_weighted_ligand_target_links：推断可能的配体和属于感兴趣基因集的基因之间的活性配体目标链接
#考虑配体的前 n 个目标与基因集之间的交叉。
active_ligand_target_links_df = best_upstream_ligands %>% 
  lapply(get_weighted_ligand_target_links,
         geneset = geneset_oi, 
         ligand_target_matrix = ligand_target_matrix, 
         n = 200) %>% 
  bind_rows() %>% 
  drop_na()
#绘制top20配体与受体的热图
active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, 
                                                                 ligand_target_matrix = ligand_target_matrix, 
                                                                 cutoff = 0.33)

order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
order_targets = active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()
rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23

vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()

p_ligand_target_network = vis_ligand_target %>% 
  make_heatmap_ggplot("Prioritized ligands","Predicted target genes", color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential")  + 
  theme(axis.text.x = element_text(face = "italic")) + 
  scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0, 0.003, 0.006))
p_ligand_target_network

#Receptors of top-ranked ligands
lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()

lr_network_top_df_large = weighted_networks_lr %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)

lr_network_top_df = lr_network_top_df_large %>% spread("from","weight",fill = 0)
lr_network_top_matrix = lr_network_top_df %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)
#按聚类顺序提取数据集中的 受配体基因矩阵
dist_receptors = dist(lr_network_top_matrix, method = "binary")
hclust_receptors = hclust(dist_receptors, method = "ward.D2")
order_receptors = hclust_receptors$labels[hclust_receptors$order]

dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
hclust_ligands = hclust(dist_ligands, method = "ward.D2")
order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]

order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix))
order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix))

# vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands]

rownames(vis_ligand_receptor_network) = order_receptors %>% make.names()
# colnames(vis_ligand_receptor_network) = order_ligands_receptor %>% make.names()
colnames(vis_ligand_receptor_network) = order_ligands %>% make.names()


p_ligand_receptor_network = vis_ligand_receptor_network %>% t() %>% make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential")
p_ligand_receptor_network

#按照更严格的筛选条件筛选受配体矩阵
#bona fide ligand-receptor interactions
lr_network_strict = lr_network %>% filter(database != "ppi_prediction_go" & database != "ppi_prediction")
ligands_bona_fide = lr_network_strict %>% pull(from) %>% unique()
receptors_bona_fide = lr_network_strict %>% pull(to) %>% unique()

lr_network_top_df_large_strict = lr_network_top_df_large %>% distinct(from,to) %>% inner_join(lr_network_strict, by = c("from","to")) %>% distinct(from,to)
lr_network_top_df_large_strict = lr_network_top_df_large_strict %>% inner_join(lr_network_top_df_large, by = c("from","to"))

lr_network_top_df_strict = lr_network_top_df_large_strict %>% spread("from","weight",fill = 0)
lr_network_top_matrix_strict = lr_network_top_df_strict %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df_strict$to)

dist_receptors = dist(lr_network_top_matrix_strict, method = "binary")
hclust_receptors = hclust(dist_receptors, method = "ward.D2")
order_receptors = hclust_receptors$labels[hclust_receptors$order]

dist_ligands = dist(lr_network_top_matrix_strict %>% t(), method = "binary")
hclust_ligands = hclust(dist_ligands, method = "ward.D2")
order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]

order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix_strict))
order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix_strict))

vis_ligand_receptor_network_strict = lr_network_top_matrix_strict[order_receptors, order_ligands_receptor]
rownames(vis_ligand_receptor_network_strict) = order_receptors %>% make.names()
colnames(vis_ligand_receptor_network_strict) = order_ligands_receptor %>% make.names()

p_ligand_receptor_network_strict = vis_ligand_receptor_network_strict %>% t() %>% make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential\n(bona fide)")
p_ligand_receptor_network_strict

#LFC in sender cells
# DE analysis for each sender cell type
# this uses a new nichenetr function - reinstall nichenetr if necessary!

get_de_table <- function(celltype, seurat_obj){
  de_table <- FindMarkers(seurat_obj, celltype, min.pct = 0.1) %>%
    rownames_to_column('gene')
  de_table <- de_table %>% as_tibble() %>% select(-p_val) %>% select(gene, avg_log2FC)
  colnames(de_table) <- c('gene', celltype)
  return(de_table)
}
#计算受配体细胞的marker基因
macrds <- readRDS("analysis/IR五个时间点/rds文件/mac_72h.rds")
Idents(macrds) <- "celltype"
#de_table_sender <- Idents(macrds) %>% levels() %>% intersect(sender_celltypes) %>% lapply(get_de_table, macrds) %>% reduce(full_join)
de_mac_marker <- read.csv("analysis/IR五个时间点/marker/mac_marker_celltype_72h.csv")
de_table_sender <- AverageExpression(macrds,features = unique(cosg_mac_marker$gene),group.by = "celltype")
de_table_sender <- data.frame(de_table_sender$RNA)
Idents(ir72) <- "celltype"
de_table_receiver <- Idents(ir72) %>% levels() %>% intersect(receiver) %>% lapply(get_de_table, ir72) %>% reduce(full_join)
Idents(ir72) <- "cellchat_mac"
#替换NA值
de_table_sender[is.na(de_table_sender)] = 0
de_table_receiver[is.na(de_table_receiver)] = 0
# Combine ligand activities with DE information
de_table_sender$gene <- rownames(de_table_sender)
ligand_activities_de = ligand_activities %>% select(test_ligand, pearson) %>% rename(ligand = test_ligand) %>% left_join(de_table_sender %>% rename(ligand = gene))
ligand_activities_de[is.na(ligand_activities_de)] = 0

# make LFC heatmap
lfc_matrix = ligand_activities_de  %>% select(-ligand, -pearson) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities_de$ligand)
rownames(lfc_matrix) = rownames(lfc_matrix) %>% make.names()

order_ligands = rev(order_ligands[order_ligands %in% rownames(lfc_matrix)])
vis_ligand_lfc = lfc_matrix[order_ligands,]
vis_ligand_lfc <- as.matrix(vis_ligand_lfc)
colnames(vis_ligand_lfc) <- sender_celltypes

# colnames(vis_ligand_lfc) = vis_ligand_lfc %>% colnames() %>% make.names()

p_ligand_lfc = vis_ligand_lfc %>% make_threecolor_heatmap_ggplot("Prioritized ligands","LFC in Sender", low_color = "midnightblue",mid_color = "white", mid = median(vis_ligand_lfc), high_color = "red",legend_position = "top", x_axis_position = "top", legend_title = "LFC") + theme(axis.text.y = element_text(face = "italic"))
p_ligand_lfc

de_table_receiver <- de_table_receiver[de_table_receiver$gene %in% order_receptors, ]
de_table_receiver[is.na(de_table_receiver)] =0
de_table_receiver <- de_table_receiver[match(order_receptors, de_table_receiver$gene), ]
de_table_receiver$gene = order_receptors
de_table_receiver[is.na(de_table_receiver)] = 0
vis_receptor_lfc <- as.data.frame(de_table_receiver)
rownames(vis_receptor_lfc) <- vis_receptor_lfc$gene
vis_receptor_lfc$gene <- NULL
vis_receptor_lfc <- as.matrix(vis_receptor_lfc)
p_receptor_lfc = vis_receptor_lfc %>% make_threecolor_heatmap_ggplot("Receptor","LFC in receiver", low_color = "midnightblue",mid_color = "white", mid = median(vis_ligand_lfc), high_color = "red",legend_position = "right", x_axis_position = "top", legend_title = "LFC") + theme(axis.text.y = element_text(face = "italic")) + coord_flip()
p_receptor_lfc

# # ligand activity heatmap
# ligand_pearson_matrix = ligand_activities %>% select(pearson) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities$test_ligand)
# 
# rownames(ligand_pearson_matrix) = rownames(ligand_pearson_matrix) %>% make.names()
# colnames(ligand_pearson_matrix) = colnames(ligand_pearson_matrix) %>% make.names()
# 
# vis_ligand_pearson = ligand_pearson_matrix[order_ligands, ] %>% as.matrix(ncol = 1) %>% magrittr::set_colnames("Pearson")
# p_ligand_pearson = vis_ligand_pearson %>% make_heatmap_ggplot("Prioritized ligands","Ligand activity", color = "darkorange",legend_position = "top", x_axis_position = "top", legend_title = "Pearson correlation coefficient\ntarget gene prediction ability)") + theme(legend.text = element_text(size = 9))
# 
# # ligand expression Seurat dotplot
# order_ligands_adapted = order_ligands
# order_ligands_adapted[order_ligands_adapted == "H2.M3"] = "H2-M3" # cf required use of make.names for heatmap visualization | this is not necessary if these ligands are not in the list of prioritized ligands!
# order_ligands_adapted[order_ligands_adapted == "H2.T23"] = "H2-T23" # cf required use of make.names for heatmap visualization | this is not necessary if these ligands are not in the list of prioritized ligands!
# rotated_dotplot = DotPlot(ir72 %>% subset(subtype %in% sender_celltypes), features = order_ligands_adapted, cols = "RdYlBu") + coord_flip() + theme(legend.text = element_text(size = 10), legend.title = element_text(size = 12)) # flip of coordinates necessary because we want to show ligands in the rows when combining all plots
# 
# figures_without_legend = cowplot::plot_grid(
#   p_ligand_pearson + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()),
#   rotated_dotplot + theme(legend.position = "none", axis.ticks = element_blank(), axis.title.x = element_text(size = 12), axis.text.y = element_text(face = "italic", size = 9), axis.text.x = element_text(size = 9,  angle = 90,hjust = 0)) + ylab("Expression in Sender") + xlab("") + scale_y_discrete(position = "right"),
#   p_ligand_lfc + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()) + ylab(""),
#   p_ligand_target_network + theme(legend.position = "none", axis.ticks = element_blank()) + ylab(""),
#   align = "hv",
#   nrow = 1,
#   rel_widths = c(ncol(vis_ligand_pearson)+6, ncol(vis_ligand_lfc) + 7, ncol(vis_ligand_lfc) + 8, ncol(vis_ligand_target)))
# 
# legends = cowplot::plot_grid(
#   ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_pearson)),
#   ggpubr::as_ggplot(ggpubr::get_legend(rotated_dotplot)),
#   ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_lfc)),
#   ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_target_network)),
#   nrow = 1,
#   align = "h", rel_widths = c(1.5, 1, 1, 1))
# 
# combined_plot = cowplot::plot_grid(figures_without_legend, legends, rel_heights = c(10,5), nrow = 2, align = "hv")
# 
# combined_plot

p4 <- plot_grid(p_ligand_lfc, p_ligand_receptor_network, NULL, p_receptor_lfc + coord_flip(), align='hv', rel_heights = c(15,10), rel_widths = c(6,15))

pdf('wala.pdf', width = 16, height=9)
print(p4)
dev.off()
