# 241012 before pre micro subset monocle
library(Seurat)
library(ggplot2)
library(monocle)
# 加载必要的库
library(devtools)
library(monocle)
library(tidyverse)
library(patchwork)
library(stringr)
library(ggplot2)
library(Seurat)
library(future)
library(rhdf5)
library(dplyr)
library(data.table)
library(Matrix)
library(rjson)
library(monocle3)
library(patchwork)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(monocle, lib.loc = "/refdir/Rlib")
# 加载Seurat对象，提取小胶质细胞亚群
# 假设scRNAsub是你已经处理过的Seurat对象
microglia_sub <- readRDS("~/r/danny-scad/rds/micro_sub.rds")
microglia_sub <- microglia_cells
# 将Seurat对象转为CellDataSet对象
data <- as(as.matrix(microglia_sub@assays$RNA$counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = microglia_sub@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
mycds <- newCellDataSet(data,
                        phenoData = pd,
                        featureData = fd,
                        expressionFamily = negbinomial.size())

# 估算大小因子和离散度
mycds <- estimateSizeFactors(mycds)
mycds <- estimateDispersions(mycds, cores=40, relative_expr = TRUE)

# 保存和加载CellDataSet对象
saveRDS(mycds, "mycds_microglia.rds")
mycds <- readRDS("mycds_microglia.rds")

# 选择发育差异表达基因（如果有特定的marker基因）
# 也可以选择clusters差异表达基因，或者Seurat选择的高变基因
diff.genes <- read.csv('cell_identify/microglia_diff_genes_wilcox.csv')
diff.genes <- subset(diff.genes, p_val_adj < 0.01)$gene
mycds <- setOrderingFilter(mycds, diff.genes)

# 使用Seurat选择的高变基因
var.genes <- VariableFeatures(microglia_sub)
mycds <- setOrderingFilter(mycds, var.genes)

# 使用Monocle选择的高变基因
disp_table <- dispersionTable(mycds)
disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
mycds <- setOrderingFilter(mycds, disp.genes)

# 比较不同基因选择方法的结果
p1 <- plot_ordering_genes(mycds)
p2 <- plot_ordering_genes(mycds)
p3 <- plot_ordering_genes(mycds, genes = disp.genes)
p1 | p2 | p3

# 创建存储拟时序结果的目录
dir.create("./pseudotime/Microglia")

# 降维并进行细胞排序
mycds <- reduceDimension(mycds, max_components = 2, method = 'DDRTree')
mycds <- orderCells(mycds)

# 保存降维和排序结果
saveRDS(mycds, "ordered_reduced_microglia.rds")

# 绘制State轨迹分布图
plot1 <- plot_cell_trajectory(mycds, color_by = "State")
ggsave("pseudotime/Microglia/State.pdf", plot = plot1, width = 6, height = 5)
ggsave("pseudotime/Microglia/State.png", plot = plot1, width = 6, height = 5)

# 绘制Cluster轨迹分布图
plot2 <- plot_cell_trajectory(mycds, color_by = "seurat_clusters")
ggsave("pseudotime/Microglia/Cluster.pdf", plot = plot2, width = 6, height = 5)
ggsave("pseudotime/Microglia/Cluster.png", plot = plot2, width = 6, height = 5)

# 绘制Celltype轨迹分布图（如果有celltype注释）
plot4 <- plot_cell_trajectory(mycds, color_by = "celltype")
ggsave("pseudotime/Microglia/Celltype.pdf", plot = plot4, width = 6, height = 5)
ggsave("pseudotime/Microglia/Celltype.png", plot = plot4, width = 6, height = 5)

# 绘制Pseudotime轨迹图
plot3 <- plot_cell_trajectory(mycds, color_by = "Pseudotime")
ggsave("pseudotime/Microglia/Pseudotime.pdf", plot = plot3, width = 6, height = 5)
ggsave("pseudotime/Microglia/Pseudotime.png", plot = plot3, width = 6, height = 5)

# 合并多个轨迹图
plotc <- plot1 | plot2 | plot3
ggsave("pseudotime/Microglia/Combination.pdf", plot = plotc, width = 10, height = 3.5)
ggsave("pseudotime/Microglia/Combination.png", plot = plotc, width = 10, height = 3.5)

# 保存细胞的拟时序信息
write.csv(pData(mycds), "pseudotime/Microglia/pseudotime.csv")

# 选择感兴趣的基因（可以是小胶质细胞相关基因）
s.genes <- c("Tmem119", "Cx3cr1", "P2ry12", "Itgax", "Apoe")

# 绘制感兴趣基因的图
p1 <- plot_genes_jitter(mycds[s.genes,], grouping = "State", color_by = "State")
p2 <- plot_genes_violin(mycds[s.genes,], grouping = "State", color_by = "State")
p3 <- plot_genes_in_pseudotime(mycds[s.genes,], color_by = "State")
plotc <- p1 | p2 | p3

# 保存基因可视化结果
ggsave("pseudotime/Microglia/plot_genes_jitter.pdf", plot = p1, width = 10, height = 3.5)
ggsave("pseudotime/Microglia/plot_genes_violin.pdf", plot = p2, width = 10, height = 3.5)
ggsave("pseudotime/Microglia/plot_genes_in_pseudotime.pdf", plot = p3, width = 10, height = 3.5)
ggsave("pseudotime/Microglia/genes_visual.png", plot = plotc, width = 8, height = 4.5)




assay <- "RNA"
group.by <- "seurat_clusters"
run_monocle3 <- function(obj,assay='RNA',ex.umap=F,group.by='seurat_clusters',label='out') {
  all <- microglia_sub
  all$ingroup <- all@meta.data[,group.by]
  expression_matrix2 <- all@assays$RNA$counts
  # 提取表达矩阵（对于Seurat v5的对象）
  expression_matrix <- all@assays[[assay]]$data
  identical(expression_matrix,expression_matrix2)
  setdiff(expression_matrix,expression_matrix2)
  
  cell_metadata <- all@meta.data
  gene_annotation <- data.frame(rownames(all), rownames(all))
  rownames(gene_annotation) <- rownames(all)
  colnames(gene_annotation) <- c("GeneSymbol", "gene_short_name")
  
  NucSeq_cds <- new_cell_data_set(expression_matrix,
                                  cell_metadata = cell_metadata,
                                  gene_metadata = gene_annotation)
  
  NucSeq_cds@int_colData$reducedDims$PCA<- all@reductions$pca@cell.embeddings
  NucSeq_cds <- reduce_dimension(NucSeq_cds, reduction_method = 'UMAP', preprocess_method = "PCA")
  
  NucSeq_cds$celltype <- NucSeq_cds$ingroup
  
  jpeg('cds.umap.jpg')
  p1 <- plot_cells(NucSeq_cds, reduction_method="UMAP", color_cells_by="celltype",show_trajectory_graph=F) + ggtitle('cds.umap')
  print(p1)
  dev.off()
  NucSeq.coembed <- all
  
  ##import  seurat umap
  if(ex.umap){
    cds.embed <- NucSeq_cds@int_colData$reducedDims$UMAP
    int.embed <- Embeddings(NucSeq.coembed, reduction = "umap")
    int.embed <- int.embed[rownames(cds.embed),]
    NucSeq_cds@int_colData$reducedDims$UMAP <- int.embed
    jpeg('int.umap.jpg')
    p1 <- plot_cells(NucSeq_cds, reduction_method="UMAP", color_cells_by="partition",show_trajectory_graph=T) + ggtitle('int.umap')
    print(p1)
    dev.off()
    
  }
  
  NucSeq_cds <- cluster_cells(NucSeq_cds, reduction_method='UMAP')
  print(length(unique(partitions(NucSeq_cds))))
  
  jpeg('partition.umap.jpg')
  p1 <- plot_cells(NucSeq_cds, show_trajectory_graph = FALSE,group_label_size = 5,color_cells_by = "partition")
  print(p1)
  dev.off()
  
  NucSeq_cds <- learn_graph(NucSeq_cds)
  
  
  jpeg('celltype.umap.jpg')
  p1 <- plot_cells(NucSeq_cds, color_cells_by = "cluster",
                   label_groups_by_cluster = FALSE, label_leaves = TRUE,
                   label_branch_points = TRUE,group_label_size = 5)
  print(p1)
  dev.off()
  
  saveRDS(NucSeq_cds,paste0(label,".rds"))
  meta <- data.frame(NucSeq_cds@colData@listData)
  write.table(meta,paste0(label,"_metadata.xls"),sep='\t',quote=F)
}
