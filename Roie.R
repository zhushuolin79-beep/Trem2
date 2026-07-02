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

scRNA <- readRDS("RDS/WAT/RNA_snn_res.0.5_QC.rds")
#Idents(scRNA) <- "celltype"
#scRNA <- subset(scRNA,ident = "mono",invert = T)
#saveRDS(scRNA,"GSE180420_肾脏2/RDS/kidney_mac_celltype2.rds")
metadata = scRNA@meta.data
head(scRNA)
summary <- table(metadata[,c('celltype','orig.ident')])

roe <- as.data.frame(ROIE(summary))

library(pheatmap)
library(RColorBrewer)
roe <- roe[rev(1:nrow(roe)), ]

#data <- roe[-6,]
pheatmap(roe, 
         scale = "none",
         display_numbers = TRUE,
         number_color = "black",
         cluster_row = FALSE,
         cluster_col = FALSE,main="Heatmap of Ro/e",
         colorRampPalette(c("#e7eb8e","#cd782d","firebrick3"))(50))


####Ro/e(邵)####
library(dplyr)
setwd("E:/数据备份/骨髓移植/")
rm(list = ls())
scRNA <- readRDS("RDS/14_sam/RNA_snn_res.0.5_QC.rds")
#细胞占比Ro/e（Ratio of Observed cell number to Expected cell number）####
#计算每个条件（sham, MI）和每个细胞类型的观察到的细胞数
observed_counts <- scRNA@meta.data %>% group_by(Sample, celltype) %>% summarise(observed = n())
observed_counts

# 计算每个条件的细胞总数
total_cells <- scRNA@meta.data %>% group_by(Sample) %>% summarise(total_cells = n())
total_cells

# 计算每种细胞类型在所有细胞中的比例
total_celltypes <- scRNA@meta.data %>% group_by(celltype) %>% summarise(total_type_count = n())
total_celltypes

# 合并细胞总数和细胞类型总数，并计算预期的细胞数
expected_counts <- left_join(observed_counts, total_cells, by = "Sample") %>%
  left_join(total_celltypes, by = "celltype") %>%
  mutate(expected = total_cells * (total_type_count / sum(total_type_count)))

expected_counts

# 计算每个条件和细胞类型的 Ro/e 比值
Ro_e_data <- expected_counts %>% mutate(Ro_e = observed / expected)
Ro_e_data

# 调整 Sample 顺序
Ro_e_data <- subset(Ro_e_data,Ro_e_data$celltype %in% un_CD45_cell)
Ro_e_data$celltype <- factor(Ro_e_data$celltype,levels = rev(un_CD45_cell))
Ro_e_data$Sample <- factor(Ro_e_data$Sample, levels = levels(scRNA@meta.data$orig.ident))

# 绘制热图并显示 Ro/e 值
ggplot(Ro_e_data, aes(x = Sample, y = celltype, fill = Ro_e)) +
  geom_tile(color = "darkgrey", linewidth =0.2) +  # 使用深灰色的格子边框
  scale_fill_gradient(low = "white", high = "darkred") +
  theme_minimal() +
  labs(x = "Sample", y = "celltype", fill = "Ro/e") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  geom_text(aes(label = round(Ro_e, 2)), color = "black", size = 3) + # 在热图中显示 Ro/e 值
  theme(axis.text.x = element_text(size = 10, face = "bold", angle = 45, hjust = 1, family = "Arial", color = "black"),   # 旋转X轴标签
        axis.text.y = element_text(size = 10, face = "bold", margin = margin(r = -2), family = "Arial", color = "black"), # Y轴字体大小，设置右边距为5
        axis.title.x = element_text(size = 12, face = "bold", vjust = -0.5, family = "Arial", color = "black"),           # X轴标题字体大小
        axis.title.y = element_text(size = 12, face = "bold", vjust = 1, family = "Arial", color = "black"),              # Y轴标题字体大小
        panel.border = element_blank(),                                                                  # 移除边框
        panel.grid = element_blank()                                                                     # 移除网格线
  )

