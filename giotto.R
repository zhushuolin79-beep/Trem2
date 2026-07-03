library(Seurat)
library(Giotto)
library(data.table)
library(dplyr)
samps <- c('IR6h')
for (name in samps){
  object_sp <- readRDS(paste0(name, '_all_info.rds'))
  spatial_locations <- as.data.table(object_sp@images$slice1@coordinates)
  colnames(spatial_locations) = c('in_tissue', 'array_row', 'array_col', 'col_pxl', 'row_pxl')
  instrs = createGiottoInstructions(python_path = '/BGFS1/home/fengk/miniconda3/envs/giotto_env/bin/python')
  raw_matrix <- object_sp@assays$Spatial@counts
  
  giotto_obj <- createGiottoObject(raw_exprs = raw_matrix,
                                   spatial_locs = spatial_locations[,.(row_pxl,-col_pxl)],
                                   instructions = instrs,
                                   cell_metadata = spatial_locations[,.(in_tissue, array_row, array_col)])
  
  DefaultAssay(object_sp) <- 'Spatial'
  giotto_obj@norm_expr <- as.matrix(object_sp@assays$Spatial@data)
  scRNAseq_markers = read.csv(file = 'ir_subtype_deg_filtered.txt', sep='\t', stringsAsFactors = F)
  top100 <- scRNAseq_markers %>% group_by(cluster) %>% top_n(100,avg_log2FC)
  
  cbind_fill <- function (rawlist) 
  {
    rown <- unique(unlist(sapply(rawlist, names)))
    matlist <- lapply(rawlist, as.matrix)
    clist <- lapply(matlist, function(x) {
      missr <- rown[!rown %in% rownames(x)]
      mat <- matrix(0, length(missr), ncol(x))
      rownames(mat) <- missr
      matc <- rbind(x, mat)
      matc <- matc[rown, ]
    })
    data.frame(clist)
  }
  
  markerlist <- list()
  for (i in unique(top100$cluster)) {
    mk <- as.vector(top100[top100$cluster == i,]$gene)
    gene <- rep(1,length(mk))
    names(gene) <- mk
    markerlist[[i]] <- gene
  }
  sig_matrix <- cbind_fill(rawlist = markerlist)
  sig_matrix <- as.matrix(sig_matrix)
  head(sig_matrix)
  
  rown <- unique(unlist(sapply(markerlist, names)))
  matlist <- lapply(markerlist, as.matrix)
  clist <- lapply(matlist, function(x) {
    missr <- rown[!rown %in% rownames(x)]
    mat <- matrix(0, length(missr), ncol(x))
    rownames(mat) <- missr
    matc <- rbind(x, mat)
    matc <- matc[rown, ]
  })
  sig_matrix <- do.call(cbind,clist)
  
  giotto_obj <- runSpatialEnrich(giotto_obj, sign_matrix = sig_matrix, enrich_method = 'PAGE')
  SPOT_frac <- data.frame(giotto_obj@spatial_enrichment$PAGE, row.names = 1, stringsAsFactors = F)
  colnames(SPOT_frac) <- gsub('-','_',colnames(SPOT_frac))
  col_order <- colnames(SPOT_frac)
  SPOT_frac <- as.data.frame(SPOT_frac[colnames(object_sp),])
  SPOT_frac <- apply(SPOT_frac, 1, function(x) {
    x <-x + abs(min(x))
    round(x/sum(x),2)
  })
  SPOT_frac <- t(SPOT_frac)
  
  SPOT_frac <- as.data.frame(SPOT_frac)
  write.table(SPOT_frac, paste0('rest_decon/', name, '_subtype_spot_frac.txt'), sep = '\t', quote = F)
}
