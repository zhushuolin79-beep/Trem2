library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)
library(patchwork)
library(tidyr)
library(stringr)
setwd("D:/课题数据/ST/")
sci_colors <- c(
  "#E41A1C", # 红色
  "#377EB8", # 蓝色
  "#4DAF4A", # 绿色
  "#984EA3", # 紫色
  "#FF7F00", # 橙色
  "#FFFF33", # 黄色
  "#A65628", # 棕色
  "#F781BF", # 粉红
  "#999999", # 灰色
  "#66C2A5", # 蓝绿色
  "#FC8D62", # 珊瑚色
  "#8DA0CB"  # 淡紫色
)
####figure1####
#空间分区
ST <- readRDS("STdata/rds/分群rds/IR0h_all_info.rds")
ST <- readRDS("STdata/rds/分群rds/IR6h_all_info.rds")
ST <- readRDS("STdata/rds/分群rds/IR24h_all_info.rds")
ST <- readRDS("STdata/rds/分群rds/IR72h_all_info.rds")
ST <- readRDS("STdata/rds/分群rds/IR30D_all_info.rds")
SpatialDimPlot(ST,label = T)
#分区功能富集
timepoint_name <- c("6h","24h","72h","30D")
zone <- c("IZ","DZ","BZ","RZ","RV")
enrich_data <- data.frame()
for(i in timepoint_name){
  allgsea <- read.csv(paste0("analysis/STdata/分区GSEA/",i,"/all_zone_GSEA.csv"))
  allgsea <- allgsea %>% arrange(match(cluster,zone))
  allgsea$cluster <- factor(allgsea$cluster,levels = zone)
  term <- read.csv(paste0("analysis/STdata/分区GSEA/",i,"/term.csv"))
  term <- term$term
  gseaplot <- subset(allgsea,allgsea$Description %in% unique(term))
  gseaplot$Description <- factor(gseaplot$Description,levels = unique(term))
  gseaplot <- gseaplot[,c(3,5,7,13)]
  gseaplot$pvalue <- -log10(gseaplot$pvalue)
  gseaplot$timepoint <- paste0(i)
  enrich_data <- rbind(enrich_data,gseaplot)
}
head(enrich_data)
# 数据预处理建议（根据实际需求调整）
enrich_data_clean <- enrich_data %>%
  mutate(
    # 保持Description原始顺序（按首次出现顺序）
    Description = factor(Description, levels = unique(Description)),
    # 保持cluster原始顺序（按首次出现顺序）
    cluster = factor(cluster, levels = unique(cluster)),
    # 保持timepoint原始顺序（按首次出现顺序）
    timepoint = factor(timepoint, levels = timepoint_name),
    # 简化长描述名称（保留前15个字符）
    Description_short = str_trunc(as.character(Description), 30, "right")
  )
head(enrich_data_clean)

# 高级分面气泡图
ggplot(enrich_data_clean, 
       aes(x = Description_short, 
           y = cluster, 
           color = enrichmentScore,
           size = pvalue)) +
  geom_point(alpha = 0.8) + # 半透明处理防止重叠
  facet_wrap(~ timepoint, nrow = 1, scales = "free_x") + # 水平排列时间点
  scale_color_gradient2(
    low = "#3182BD",   # 蓝色表示负值
    mid = "#F0F0F0",   # 灰色中间值
    high = "#DE2D26",  # 红色表示正值
    midpoint = 0,      # 设置中性点
    limits = c(-1, 1)  # 根据实际数据范围调整
  ) +
  scale_size_continuous(
    range = c(1, 6),   # 点大小范围
    breaks = c(1, 5, 10), # 显著性断点设置
    labels = c("1", "1e-5", "1e-10") # 显示原始p值级别
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1),
    panel.grid.major = element_line(color = "grey90"),
    strip.background = element_rect(fill = "#A65628"),
    strip.text = element_text(color = "white", face = "bold"),
    legend.position = "bottom",
    legend.box = "horizontal"
  ) +
  labs(
    x = "Biological Pathways",
    y = "Spatial Region",
    color = "Enrichment Score",
    size = "-log10(p)"
  ) +
  guides(
    color = guide_colorbar(barwidth = 15, barheight = 1),
    size = guide_legend(nrow = 1)
  )
#修改
library(dplyr)
library(ggplot2)
library(stringr)

timepoint_name <- c("6h","24h","72h","30D")
zone <- c("IZ","DZ","BZ","RZ","RV")
enrich_data <- data.frame()

for(i in timepoint_name){
  allgsea <- read.csv(paste0("analysis/STdata/分区GSEA/",i,"/all_zone_GSEA.csv"))
  allgsea <- allgsea %>% arrange(match(cluster,zone))
  allgsea$cluster <- factor(allgsea$cluster,levels = zone)
  
  term <- read.csv(paste0("analysis/STdata/分区GSEA/",i,"/term.csv"))
  term_list <- unique(term$term)
  
  gseaplot <- subset(allgsea, allgsea$Description %in% term_list)
  gseaplot$Description <- factor(gseaplot$Description, levels = term_list)
  
  # 推荐修改：使用列名直接提取，确保准确抓取 NES 和 setSize
  gseaplot <- gseaplot %>% 
    select(cluster, Description, NES, setSize, p.adjust) %>%
    # 如果你也想在这里过滤掉不显著的通路，可以取消下一行的注释
    # filter(p.adjust < 0.05) %>% 
    mutate(timepoint = paste0(i))
  
  enrich_data <- rbind(enrich_data, gseaplot)
}

head(enrich_data)

# 数据预处理（保持不变，非常好）
enrich_data_clean <- enrich_data %>%
  mutate(
    Description = factor(Description, levels = unique(Description)),
    cluster = factor(cluster, levels = unique(cluster)),
    timepoint = factor(timepoint, levels = timepoint_name),
    Description_short = str_trunc(as.character(Description), 30, "right")
  )
head(enrich_data_clean)

# 高级分面气泡图
ggplot(enrich_data_clean, 
       aes(x = Description_short, 
           y = cluster, 
           color = NES,         # 修改 1：颜色映射为 NES
           size = setSize)) +   # 修改 2：大小映射为 setSize (基因集大小)
  geom_point(alpha = 0.8) + 
  facet_wrap(~ timepoint, nrow = 1, scales = "free_x") + 
  scale_color_gradient2(
    low = "#3182BD",   
    mid = "#F0F0F0",   
    high = "#DE2D26",  
    midpoint = 0
    # 注意：已删除 limits = c(-1, 1)，以防真实的 NES (如 1.8) 无法显示颜色
  ) +
  scale_size_continuous(
    range = c(2, 6),
    breaks = c(50,100,200)
    # 调整点大小的视觉范围，让大基因集和小基因集对比明显
    # 已删除原先用于 pvalue 的 breaks 和 labels，让系统根据真实的基因数自动生成图例
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.grid.major = element_line(color = "grey90"),
    strip.background = element_rect(fill = "#A65628"),
    strip.text = element_text(color = "white", face = "bold"),
    legend.position = "bottom",
    legend.box = "horizontal"
  ) +
  labs(
    x = "Biological Pathways",
    y = "Spatial Region",
    color = "NES",            # 修改 3：图例标签更新为 NES
    size = "Gene Set Size"    # 修改 4：图例标签更新为基因集大小
  ) +
  guides(
    color = guide_colorbar(barwidth = 15, barheight = 1),
    size = guide_legend(nrow = 1)
  )
#空间celltype映射结果分区统计
library(ggradar)
data <- read.csv("analysis/STdata/反卷积/ecsub/dotplot_ECsub_repair_30D - 2.csv")
#colnames(data) <- c(0,0:6)
p5 <- ggradar(
  data,             # 分组颜色
  grid.min = 0,                    # 最小值
  grid.mid = 0.025,                  # 中间值
  grid.max = 0.05,                    # 最大值
  legend.position = "bottom",      # 图例位置
  axis.label.size = 4,             # 轴标签字体大小
  grid.label.size = 0,             # 网格标签字体大小
  background.circle.colour = "gray",# 背景圆圈颜色
  background.circle.transparency = 0.1, # 背景圆圈透明度
  group.colours = sci_colors[1:7],
  group.line.width = 1,          # 分组线宽度
  group.point.size = 2,            # 分组点大小
  fill = T,                     # 是否填充
  fill.alpha = 0.1                 # 填充透明度
)+
  annotate("text", x = 0, y = 0, label = "0", size = 5, color = "black") +  # 最小值标签
  annotate("text", x = 0, y = 0.025, label = "0.025", size = 5, color = "black") +  # 中间值标签
  annotate("text", x = 0, y = 0.05, label = "0.05", size = 5, color = "black")# 最大值标签
p1
# 横向排列合并（5图一行）
combined_plot <- p1 + p2 + p3 + p4 + p5 + 
  plot_layout(nrow = 1)  # 关键合并参数

#亚群分析
scRNA <- readRDS("analysis/IR五个时间点/rds文件/MAC_select_combined.rds")
DefaultAssay(scRNA) <- "RNA"
VlnPlot(scRNA,features = "Trem2",pt.size = 0)
DimPlot(scRNA,label = T,pt.size = 1)
metadata=scRNA@meta.data
summary <- table(metadata[,c('celltype2','orig.ident')])

ROIE <- function(crosstab){
  ## Calculate the Ro/e value from the given crosstab
  ##
  ## Args:
  #' @crosstab: the contingency table of given distribution
  ##
  ## Return:
  ## The Ro/e matrix 
  rowsum.matrix <- matrix(0, nrow = nrow(crosstab), ncol = ncol(crosstab))
  rowsum.matrix[,1] <- rowSums(crosstab)
  colsum.matrix <- matrix(0, nrow = ncol(crosstab), ncol = ncol(crosstab))
  colsum.matrix[1,] <- colSums(crosstab)
  allsum <- sum(crosstab)
  roie <- divMatrix(crosstab, rowsum.matrix %*% colsum.matrix / allsum)
  row.names(roie) <- row.names(crosstab)
  colnames(roie) <- colnames(crosstab)
  return(roie)
}
divMatrix <- function(m1, m2){
  ## Divide each element in turn in two same dimension matrixes
  ##
  ## Args:
  #' @m1: the first matrix
  #' @m2: the second matrix
  ##
  ## Returns:
  ## a matrix with the same dimension, row names and column names as m1. 
  ## result[i,j] = m1[i,j] / m2[i,j]
  dim_m1 <- dim(m1)
  dim_m2 <- dim(m2)
  if( sum(dim_m1 == dim_m2) == 2 ){
    div.result <- matrix( rep(0,dim_m1[1] * dim_m1[2]) , nrow = dim_m1[1] )
    row.names(div.result) <- row.names(m1)
    colnames(div.result) <- colnames(m1)
    for(i in 1:dim_m1[1]){
      for(j in 1:dim_m1[2]){
        div.result[i,j] <- m1[i,j] / m2[i,j]
      }
    }   
    return(div.result)
  }
  else{
    warning("The dimensions of m1 and m2 are different")
  }
}

roe <- as.data.frame(ROIE(summary))

library(pheatmap)
library(RColorBrewer)
# 生成显著性标记矩阵
significance_matrix <- matrix(
  ifelse(roe > 2, "**",
         ifelse(roe > 1, "*", "")),
  nrow = nrow(roe),
  dimnames = dimnames(roe)
)
#data <- roe[-6,]
pheatmap(roe, 
         scale = "row",
         display_numbers = significance_matrix,
         fontsize_number = 20,
         number_color = "black",
         cluster_row = FALSE,
         cluster_col = FALSE,main="Heatmap of Ro/e",
         colorRampPalette(c("#3182BD","#F0F0F0","#DE2D26"))(50))
#亚群功能
library(ggplot2)
library(dplyr)
library(tidyr)

# 加载必要的包
library(ggplot2)
library(dplyr)
library(tidyr)

# 1. 读取数据
allgsea <- read.csv("analysis/IR五个时间点/EC GSEA/all_ECsub_GSEA.csv")
term <- read.csv("analysis/IR五个时间点/EC GSEA/term.csv")
term_list <- term$term.top1.

# 2. 数据过滤与转换
gseaplot <- allgsea %>%
  filter(Description %in% term_list) %>%
  dplyr::select(cluster, Description, NES, setSize, p.adjust)

# 3. 修复 Y 轴顺序 (严格绑定至 term.csv 的顺序，rev() 确保第一行在图表最上方)
gseaplot$Description <- factor(gseaplot$Description, levels = unique(term_list))

# 4. 绘制气泡图
ggplot(gseaplot, aes(x = cluster, y = Description)) +
  # 气泡大小映射为基因数，颜色映射回 NES
  geom_point(aes(size = setSize, color = NES)) +
  # 使用双向渐变色：蓝(负) - 白(0) - 红(正)，完美契合 GSEA 的上下调含义
  scale_color_gradient2(low = "#3182BD", mid = "#F0F0F0", high = "#DE2D26", midpoint = 0) +
  scale_size_continuous(range = c(2, 8)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title = element_blank(),
    panel.grid.major = element_line(color = "grey90", linetype = "dashed"),
    panel.grid.minor = element_blank()
  ) +
  labs(
    size = "Gene Set Size",
    color = "NES"
  )
ggsave("GSEA_Dotplot_EC.pdf", width = 8, height = 6)

#空间细胞亚型映射雷达图
data <- read.csv("analysis/STdata/反卷积/ecsub/dotplot_ECsub_repair_30D - 2.csv")
#colnames(data) <- c("cluster",0:8)
p5 <- ggradar(
  data,             # 分组颜色
  grid.min = 0,                    # 最小值
  grid.mid = 0.025,                  # 中间值
  grid.max = 0.05,                    # 最大值
  legend.position = "bottom",      # 图例位置
  axis.label.size = 4,             # 轴标签字体大小
  grid.label.size = 0,             # 网格标签字体大小
  background.circle.colour = "gray",# 背景圆圈颜色
  background.circle.transparency = 0.1, # 背景圆圈透明度
  group.colours = sci_colors[1:7],
  group.line.width = 1.5,          # 分组线宽度
  group.point.size = 3,            # 分组点大小
  fill = T,                     # 是否填充
  fill.alpha = 0.1                 # 填充透明度
)+
  annotate("text", x = 0, y = 0, label = "0", size = 5, color = "black") +  # 最小值标签
  annotate("text", x = 0, y = 0.025, label = "0.025", size = 5, color = "black") +  # 中间值标签
  annotate("text", x = 0, y = 0.05, label = "0.05", size = 5, color = "black")# 最大值标签
# 横向排列合并（5图一行）
combined_plot <- p1 + p2 + p3 + p4 + p5 + 
  plot_layout(nrow = 1)  # 关键合并参数
combined_plot
p1
##MISTy分析
#柱状图统计
data <- read.csv("analysis/STdata/misty/修复细胞类型5-3/misty_all_result.csv")
head(data)
# 计算每个Predictor在不同zone的importance总和
result <- data %>%
  filter(view == "intra_5", group == "72h") %>%  # 确保过滤条件正确
  group_by(Predictor, zone) %>%
  summarise(total_importance = sum(Importance,na.rm = TRUE), .groups = "drop")
result
# 可视化（可选）
library(ggplot2)
ggplot(result, aes(x = reorder(Predictor, -total_importance), y = total_importance)) +
  geom_col(fill = "steelblue") +
  labs(title = "Total Importance by Predictor (6H)",
       x = "Predictor", y = "Total Importance") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

predictor_name <- c("Trem2.mac","Lyve1.mac","MYO","Ki67.CF","Sca1hi.CF","EC.cap","EndoMT","EC.proliferation")
result$zone <- factor(result$zone,levels = c("IZ","DZ","BZ","RZ","RV"))
result$Predictor <- factor(result$Predictor,levels = predictor_name)
ggplot(result, aes(x = Predictor, y = total_importance, fill = zone)) +
  geom_col(position = "dodge") +  # 并排显示不同 zone
  facet_wrap(~zone, scales = "free_x") +  # 按 zone 分面
  labs(
    title = "MISTy--Total Importance by Predictor and Zone (Juxta_view, 72H)",
    x = "Predictor",
    y = "Total Importance"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  NoLegend()
####figure1-sup1####
#dimplot and dotplot  
ST <- readRDS("analysis/STdata/rds/分群rds/IR30D_all_info.rds")
marker <- read.csv("analysis/STdata/marker/ST_30D.csv")
top3 <- marker %>% group_by(cluster) %>% top_n(3,avg_log2FC)
DimPlot(ST,label = T,pt.size = 1) + NoLegend()
DotPlot(ST ,features = unique(top3$gene),cols = c("#3182BD","#DE2D26"),col.min = -1,col.max = 2)+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))
#celltype mapping
ST <- readRDS("analysis/STdata/rds/映射rds/IR24H_with_major_frac.rds")
DefaultAssay(ST) <- 'major_frac'
SpatialFeaturePlot(ST,features = c("EC"),min.cutoff = 0,max.cutoff = 0.2)
#scRNA-seq celltype
scRNA <- readRDS("analysis/IR五个时间点/rds文件/IR0h_30d_修改.rds")
VlnPlot(scRNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), 
        group.by = "orig.ident",ncol = 3,pt.size = 0)
FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",
               group.by = "orig.ident",pt.size = 1)
DimPlot(scRNA, reduction = "pca", group.by = "orig.ident",pt.size = 1)
DimPlot(scRNA, reduction = "umap", group.by = "celltype",label = T,pt.size = 1)
DimPlot(scRNA, reduction = "umap", group.by = "celltype",
        split.by = "orig.ident",label = F,pt.size = 1)
#Ro/e
metadata=scRNA@meta.data
summary <- table(metadata[,c('celltype','orig.ident')])

ROIE <- function(crosstab){
  ## Calculate the Ro/e value from the given crosstab
  ##
  ## Args:
  #' @crosstab: the contingency table of given distribution
  ##
  ## Return:
  ## The Ro/e matrix 
  rowsum.matrix <- matrix(0, nrow = nrow(crosstab), ncol = ncol(crosstab))
  rowsum.matrix[,1] <- rowSums(crosstab)
  colsum.matrix <- matrix(0, nrow = ncol(crosstab), ncol = ncol(crosstab))
  colsum.matrix[1,] <- colSums(crosstab)
  allsum <- sum(crosstab)
  roie <- divMatrix(crosstab, rowsum.matrix %*% colsum.matrix / allsum)
  row.names(roie) <- row.names(crosstab)
  colnames(roie) <- colnames(crosstab)
  return(roie)
}
divMatrix <- function(m1, m2){
  ## Divide each element in turn in two same dimension matrixes
  ##
  ## Args:
  #' @m1: the first matrix
  #' @m2: the second matrix
  ##
  ## Returns:
  ## a matrix with the same dimension, row names and column names as m1. 
  ## result[i,j] = m1[i,j] / m2[i,j]
  dim_m1 <- dim(m1)
  dim_m2 <- dim(m2)
  if( sum(dim_m1 == dim_m2) == 2 ){
    div.result <- matrix( rep(0,dim_m1[1] * dim_m1[2]) , nrow = dim_m1[1] )
    row.names(div.result) <- row.names(m1)
    colnames(div.result) <- colnames(m1)
    for(i in 1:dim_m1[1]){
      for(j in 1:dim_m1[2]){
        div.result[i,j] <- m1[i,j] / m2[i,j]
      }
    }   
    return(div.result)
  }
  else{
    warning("The dimensions of m1 and m2 are different")
  }
}

roe <- as.data.frame(ROIE(summary))

library(pheatmap)
library(RColorBrewer)
# 生成显著性标记矩阵
significance_matrix <- matrix(
  ifelse(roe > 2, "**",
         ifelse(roe > 1, "*", "")),
  nrow = nrow(roe),
  dimnames = dimnames(roe)
)
#data <- roe[-6,]
pheatmap(roe, 
         scale = "row",
         display_numbers = significance_matrix,
         fontsize_number = 20,
         number_color = "black",
         cluster_row = FALSE,
         cluster_col = FALSE,main="Heatmap of Ro/e",
         colorRampPalette(c("#3182BD","#F0F0F0","#DE2D26"))(50))
#vlnplot
VlnPlot(scRNA,features = c("Cd68","Pdgfra","S100a9","Cd79a","H2-Ab1","Cd3d",
                           "Pecam1","Rgs5","Acta2","Myl2","Cd79b"),pt.size = 0,stack = T) + NoLegend()
####Figure1-sup2####
#celltype-sub
scRNA <- readRDS("analysis/IR五个时间点/rds文件/MAC_select_combined.rds")
DimPlot(scRNA, reduction = "umap", group.by = "",
        split.by = "orig.ident",label = F,pt.size = 1)
marker <- read.csv("analysis/IR五个时间点/marker/mac marker/mac_marker_celltype.csv")
top3 <- marker %>% group_by(cluster) %>% top_n(3,avg_log2FC)
genename <- top3$gene
genename[c(15,17,23)] <- c("Ltb4r1","Ifit2","Mki67")
DotPlot(scRNA ,features = genename,cols = c("#3182BD","#DE2D26"),col.min = -1,col.max = 2)+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))
#ST-celltypesub mapping
repair_celltype <- c("Trem2-mac","Lyve1-mac","MYO","Sca1hi-CF","Ki67-CF","EC.cap","EndoMT",'EC.proliferation')
ST <- readRDS("analysis/STdata/rds/映射rds/IR24h_with_major_frac.rds")
DefaultAssay(ST) <- 'frac'
SpatialFeaturePlot(ST,features = repair_celltype,min.cutoff = 0,max.cutoff = 0.1,ncol = 8)
SpatialFeaturePlot(ST,features = repair_celltype[1],min.cutoff = 0,max.cutoff = 0.1)

####Figure2####
setwd("D:/课题数据/ST/")
scRNA <- readRDS("analysis/公共数据/GSE210159_心脏2/RDS/heart_mac_celltype.rds")
DimPlot(scRNA,label = T,pt.size = 1)
VlnPlot(scRNA,"Dhcr7")
#vlnplot
genelist <- c("Trem2","H2-Eb1","Lyve1","Ltb4r1","Ifit2","Mki67")
VlnPlot(scRNA,features = genelist,pt.size = 0,stack = T)
#Ro/e
metadata=scRNA@meta.data
head(scRNA)
summary <- table(metadata[,c('celltype','group')])

ROIE <- function(crosstab){
  ## Calculate the Ro/e value from the given crosstab
  ##
  ## Args:
  #' @crosstab: the contingency table of given distribution
  ##
  ## Return:
  ## The Ro/e matrix 
  rowsum.matrix <- matrix(0, nrow = nrow(crosstab), ncol = ncol(crosstab))
  rowsum.matrix[,1] <- rowSums(crosstab)
  colsum.matrix <- matrix(0, nrow = ncol(crosstab), ncol = ncol(crosstab))
  colsum.matrix[1,] <- colSums(crosstab)
  allsum <- sum(crosstab)
  roie <- divMatrix(crosstab, rowsum.matrix %*% colsum.matrix / allsum)
  row.names(roie) <- row.names(crosstab)
  colnames(roie) <- colnames(crosstab)
  return(roie)
}
divMatrix <- function(m1, m2){
  ## Divide each element in turn in two same dimension matrixes
  ##
  ## Args:
  #' @m1: the first matrix
  #' @m2: the second matrix
  ##
  ## Returns:
  ## a matrix with the same dimension, row names and column names as m1. 
  ## result[i,j] = m1[i,j] / m2[i,j]
  dim_m1 <- dim(m1)
  dim_m2 <- dim(m2)
  if( sum(dim_m1 == dim_m2) == 2 ){
    div.result <- matrix( rep(0,dim_m1[1] * dim_m1[2]) , nrow = dim_m1[1] )
    row.names(div.result) <- row.names(m1)
    colnames(div.result) <- colnames(m1)
    for(i in 1:dim_m1[1]){
      for(j in 1:dim_m1[2]){
        div.result[i,j] <- m1[i,j] / m2[i,j]
      }
    }   
    return(div.result)
  }
  else{
    warning("The dimensions of m1 and m2 are different")
  }
}

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
#相关性分析
library(ggpmisc)
library(ggpubr)
library(dplyr)
refdata <- readRDS("D:/课题数据/ST/analysis/IR五个时间点/rds文件/MAC_select_combined.rds")
marker <- read.csv("D:/课题数据/ST/analysis/IR五个时间点/marker/mac marker/mac_marker_celltype_all.csv")
marker <- subset(marker, marker$avg_log2FC > 0 & marker$p_val_adj < 0.05 & marker$cluster == "Trem2_mac")
gene <- marker$gene
gene_exp1 <- data.frame(AverageExpression(refdata,features = gene,group.by = "celltype",assays = "RNA"))
colnames(gene_exp1) <- paste0("Heart_IR_",colnames(gene_exp1))
#gene_exp1$gene <- rownames(gene_exp1)
head(gene_exp1)

scRNA <- readRDS("D:/课题数据/ST/analysis/公共数据/GSE227075_下肢/RDS文件/climb_mac_celltype.rds")
gene_exp2 <- data.frame(AverageExpression(scRNA,features = gene,group.by = "celltype",assays = "RNA"))
colnames(gene_exp2) <- paste0("climb_IR_",colnames(gene_exp2))
#gene_exp2$gene <- rownames(gene_exp2)
head(gene_exp2)

gene_exp1 <- subset(gene_exp1,rownames(gene_exp1) %in% rownames(gene_exp2))
data_plot <- cbind(gene_exp1,gene_exp2)
cor=round(cor(data_plot$Heart_IR_RNA.Trem2_mac,
              data_plot$climb_IR_RNA.Trem2_mac),
          2)
ggplot(log1p(data_plot), aes(x = Heart_IR_RNA.Trem2_mac, y = climb_IR_RNA.Trem2_mac)) +
  geom_point(color="#6baed6")+
  geom_smooth(method = "lm", formula = y~x, color = "#756bb1", fill = "#cbc9e2")+ 
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  stat_cor(size = 8)
#功能富集
go_data <- read.csv("analysis/公共数据/GSE227075_下肢/GSEA分析72H/Trem2_mac_GSEA.csv")
term <- read.csv("analysis/公共数据/GSE227075_下肢/GSEA分析72H/term.csv")
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
####figure2-sup2####
setwd("D:/课题数据/ST/")
scRNA <- readRDS("analysis/公共数据/GSE227075_下肢/RDS文件/climb_celltype.rds")
DimPlot(scRNA,label = T,pt.size = 1,group.by = "celltype")
#亚群Roe
scRNA <- readRDS("analysis/公共数据/GSE227075_下肢/RDS文件/climb_mac_celltype.rds")
ROIE <- function(crosstab){
  ## Calculate the Ro/e value from the given crosstab
  ##
  ## Args:
  #' @crosstab: the contingency table of given distribution
  ##
  ## Return:
  ## The Ro/e matrix 
  rowsum.matrix <- matrix(0, nrow = nrow(crosstab), ncol = ncol(crosstab))
  rowsum.matrix[,1] <- rowSums(crosstab)
  colsum.matrix <- matrix(0, nrow = ncol(crosstab), ncol = ncol(crosstab))
  colsum.matrix[1,] <- colSums(crosstab)
  allsum <- sum(crosstab)
  roie <- divMatrix(crosstab, rowsum.matrix %*% colsum.matrix / allsum)
  row.names(roie) <- row.names(crosstab)
  colnames(roie) <- colnames(crosstab)
  return(roie)
}
divMatrix <- function(m1, m2){
  ## Divide each element in turn in two same dimension matrixes
  ##
  ## Args:
  #' @m1: the first matrix
  #' @m2: the second matrix
  ##
  ## Returns:
  ## a matrix with the same dimension, row names and column names as m1. 
  ## result[i,j] = m1[i,j] / m2[i,j]
  dim_m1 <- dim(m1)
  dim_m2 <- dim(m2)
  if( sum(dim_m1 == dim_m2) == 2 ){
    div.result <- matrix( rep(0,dim_m1[1] * dim_m1[2]) , nrow = dim_m1[1] )
    row.names(div.result) <- row.names(m1)
    colnames(div.result) <- colnames(m1)
    for(i in 1:dim_m1[1]){
      for(j in 1:dim_m1[2]){
        div.result[i,j] <- m1[i,j] / m2[i,j]
      }
    }   
    return(div.result)
  }
  else{
    warning("The dimensions of m1 and m2 are different")
  }
}

metadata=scRNA@meta.data
head(scRNA)
summary <- table(metadata[,c('celltype','group')])
roe <- as.data.frame(ROIE(summary))

library(pheatmap)
library(RColorBrewer)
# 生成显著性标记矩阵
significance_matrix <- matrix(
  ifelse(roe > 2, "**",
         ifelse(roe > 1, "*", "")),
  nrow = nrow(roe),
  dimnames = dimnames(roe)
)
#data <- roe[-6,]
pheatmap(roe, 
         scale = "row",
         display_numbers = significance_matrix,
         fontsize_number = 20,
         number_color = "black",
         cluster_row = FALSE,
         cluster_col = FALSE,main="Heatmap of Ro/e",
         colorRampPalette(c("#3182BD","#F0F0F0","#DE2D26"))(50))
#巨噬细胞各亚群功能
library(tibble)
allgsea <- read.csv("analysis/公共数据/GSE210159_心脏2/GSEA分析/all_heart_mac_GSEA.csv")
allgsea <-  allgsea %>% arrange(desc(allgsea$enrichmentScore))
term <- read.csv("analysis/公共数据/GSE210159_心脏2/GSEA分析/term.csv")
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
         main = "GSEA of macsub in heart_MI",
         colorRampPalette(c("#3182BD","#F0F0F0","#DE2D26"))(50))

####figure4####
#空间通讯网络
IR0H <- readRDS("analysis/STdata/空间cellchat/IR0h/cellchat.rds")
IR6H <- readRDS("analysis/STdata/空间cellchat/IR6h/cellchat.rds")
IR24H <- readRDS("analysis/STdata/空间cellchat/IR24h/cellchat.rds")
IR72H <- readRDS("analysis/STdata/空间cellchat/IR72h/cellchat.rds")
IR30D <- readRDS("analysis/STdata/空间cellchat/IR30D/cellchat.rds")
cellchat <- IR30D
groupSize <- as.numeric(table(cellchat@idents))
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, title.name = "Number of interactions",arrow.size = 0.5)
#组合分析
object.list <- list(IR6H = IR6H,
                    IR24H = IR24H,
                    IR72H = IR72H,
                    IR30D = IR30D)
cellchat <- mergeCellChat(object.list, add.names = names(object.list),cell.prefix = TRUE)
groupSize <- as.numeric(table(cellchat@idents))
cellchat@idents <- factor(cellchat@idents,levels = celltype_all )
compareInteractions(cellchat, show.legend = F, group = c(1,2,3,4,5))

netVisual_bubble(cellchat, 
                 signaling = "NOTCH",
                 sources.use = c(1), 
                 targets.use = c(2,3,4,5), 
                 comparison = c(1,2,3,4))
# 1. 运行绘图代码，去除 source 和 target 限制以获取 NOTCH 通路全集，并将结果赋值给变量 p
p <- netVisual_bubble(cellchat, 
                      signaling = "NOTCH",
                      # 不设置 sources.use 和 targets.use，代表提取所有细胞群之间的互作
                      comparison = c(1,2,3,4))

# 2. 从 ggplot 对象中提取底层数据框
notch_all_data <- p$data

# 3. 查看数据前几行
head(notch_all_data)

# 4. 导出为 CSV 表格，方便查看和放入文章附件
write.csv(notch_all_data, file = "NOTCH_All_Interactions_Merged.csv", row.names = FALSE)
#单细胞通讯韦恩图
library(readxl)
library(VennDiagram)
CF <- read_xlsx("analysis/IR五个时间点/三种细胞通讯取交集/72H/cross_72H_LR_pair_ CF.xlsx")
SMC <- read_xlsx("analysis/IR五个时间点/三种细胞通讯取交集/72H/cross_72H_LR_pair_ SMC.xlsx")
Pericyte <- read_xlsx("analysis/IR五个时间点/三种细胞通讯取交集/72H/cross_72H_LR_pair_ pericyte.xlsx")
EC <- read_xlsx("analysis/IR五个时间点/三种细胞通讯取交集/72H/cross_72H_LR_pair_ EC.xlsx")
CM <- read_xlsx("analysis/IR五个时间点/三种细胞通讯取交集/72H/cross_72H_LR_pair_ CM.xlsx")
cellphone = unique(na.omit(c(CF$`cellphone_list[[i]]`,
                             SMC$`cellphone_list[[i]]`,
                             Pericyte$`cellphone_list[[i]]`,
                             EC$`cellphone_list[[i]]`,
                             CM$`cellphone_list[[i]]`)))
cellchat = unique(na.omit(c(CF$`cellchat_list[[i]]`,
                            SMC$`cellchat_list[[i]]`,
                            Pericyte$`cellchat_list[[i]]`,
                            EC$`cellchat_list[[i]]`,
                            CM$`cellchat_list[[i]]`)))
nichenet = unique(na.omit(CF$Nichenet))
cross = unique(na.omit(SMC$交集))
venn.diagram(x=list(cellphone,cellchat,nichenet),
             
             scaled = T, # 根据比例显示大小
             
             alpha= 0.4, #透明度
             
             lwd=1,lty=1,col=sci_colors[1:3], #圆圈线条粗细、形状、颜色；1 实线, 2 虚线, blank无线条
             
             label.col ='black' , # 数字颜色abel.col=c('#FFFFCC','#CCFFFF',......)根据不同颜色显示数值颜色
             
             cex = 2, # 数字大小
             
             fontface = "bold",  # 字体粗细；加粗bold
             
             fill=sci_colors[1:3], # 填充色 配色https://www.58pic.com/
             
             category.names = c("cellphone","cellchat","nichenet") , #标签名
             
             cat.dist = 0.02, # 标签距离圆圈的远近
             
             cat.pos = c(-120, -240, -180), # 标签相对于圆圈的角度cat.pos = c(-10, 10, 135)
             
             cat.cex = 2, #标签字体大小
             
             cat.fontface = "bold",  # 标签字体加粗
             
             cat.col='black' ,   #cat.col=c('#FFFFCC','#CCFFFF',.....)根据相应颜色改变标签颜色
             
             cat.default.pos = "outer",  # 标签位置, outer内;text 外
             
             output=T,
             
             filename="figure/figure4/细胞通讯分析/venn.pdf",# 文件保存
             
             resolution = 200,  # 分辨率
             
             #compression = "lzw"# 压缩算法
             
)
#stLearn
library(pheatmap)
ST_6h_DZ <- read.csv("analysis/STdata/stlearn/重新统计/IR72h_IZ.csv",row.names = 1)
bk <- c(seq(0,50,by=10))
pheatmap(ST_6h_DZ,
         show_rownames = T,
         show_colnames = T,
         cluster_rows = F,
         cluster_cols = F,
         scale = "none",
         color = c(colorRampPalette(colors = c("white","red"))(length(bk))),
         legend_breaks=seq(0,50,10),
         breaks=bk)
#scRNA-seq
scRNA <- readRDS("analysis/IR五个时间点/rds文件/MAC_select_combined.rds")
DefaultAssay(scRNA) <- "RNA"
Idents(scRNA) <- "orig.ident"
scRNAsub <- subset(scRNA,ident = "IR72H")
Idents(scRNAsub) <- "celltype"
VlnPlot(scRNAsub,c("Jag1","Spp1","Tnfsf12","Igf1"),pt.size = 0.1,ncol = 2)
VlnPlot(scRNAsub,c("Notch1","Notch2","Itgav","Tnfrsf12a","Igf1r"),pt.size = 0.1,ncol = 2)
FeaturePlot(scRNAsub,"Notch1",label = T,min.cutoff = 1,pt.size = 1,order = T)
VlnPlot(scRNAsub,"Jag1",pt.size = 0.1)
#ST
ST72H <- readRDS("analysis/STdata/rds/分群rds/IR72h_all_info.rds")
SpatialFeaturePlot(ST72H,features = c("Jag1","Notch1"),min.cutoff = 0,max.cutoff = 1)
SpatialFeaturePlot(ST72H,features = c("Spp1","Tnfsf12","Igf1"),min.cutoff = 0,max.cutoff = 1)
SpatialFeaturePlot(ST72H,features = c("Notch2","Itgav","Tnfrsf12a","Igf1r"),
                   min.cutoff = 0,max.cutoff = 1)

####figure4-sup1####
#cellchat
data <- read.csv("analysis/IR五个时间点/cellchat/巨噬细胞各亚群与其他细胞类型/ir72h/net.csv")
head(data)
data <- subset(data,data$source == "Trem2_mac" & data$target %in% c("CM","CF","SMC","EC","Pericyte"))
data <- data[,c(4,5,6)]
data_pheat <- data[!duplicated(data[1:2]),]
data_pheatmap <- spread(data_pheat,key = "receptor",value = "prob")
#key为需要转置的列名，value为需要填充的值
rownames(data_pheatmap) <- data_pheatmap$ligand
data_pheatmap <- data_pheatmap[,-1]
data_pheatmap[is.na(data_pheatmap)] <- 0
data_pheatmap[data_pheatmap >= 0.002] = 0.01
p_data_pheatmap = as.matrix(data_pheatmap) %>% 
  make_heatmap_ggplot("Ligands","Receptors", 
                      color = "mediumvioletred", 
                      x_axis_position = "bottom",
                      legend_title = "L-R pair",
                      legend_position = "right")
p_data_pheatmap
#nichenet
data <- read.csv("analysis/IR五个时间点/nichenet分析/nichenet分析结果/72h/L_R.csv")
rownames(data) <- data$X
data <- data[,-1]
data_pheatmap = as.matrix(data) %>% 
  make_heatmap_ggplot("Ligands","Receptors", 
                      color = "mediumvioletred", 
                      x_axis_position = "bottom",
                      legend_title = "L-R pair",
                      legend_position = "right")
data_pheatmap
#cellphoneDB
data <- read.csv("analysis/IR五个时间点/cellphoneDB/分析结果/cellphonedb_result/72h_cout/Trem2_unimmu.csv")
head(data)
rownames(data) <- data$interacting_pair
data <- data[,-1]
miss <- c()
for(i in 1:nrow(data)) {
  if(length(which(is.na(data[i,]))) > 0.7*ncol(data)) 
    miss <- append(miss,i) 
}
data_pheatmap <- data[-miss,]

data_pheatmap[is.na(data_pheatmap)] <- 0
data_pheatmap = as.matrix(data_pheatmap) %>% 
  make_heatmap_ggplot("Ligands","Receptors", 
                      color = "mediumvioletred", 
                      x_axis_position = "bottom",
                      legend_title = "L-R pair",
                      legend_position = "right")
data_pheatmap
#空间共定位
IR72H <- readRDS("analysis/STdata/空间cellchat/IR72h/cellchat.rds")
pair <- c("SPP1_ITGAV_ITGB1","TNFSF12_TNFRSF12A","IGF1_IGF1R","JAG1_NOTCH1","JAG1_NOTCH2")
spatialFeaturePlot(IR72H, pairLR.use = pair, 
                    do.binary = T, 
                   cutoff = 0.05, enriched.only = F, 
                   color.heatmap = "Reds", direction = 1,
                   ncol = 5)
