rm(list = ls())
getwd()
scRNA <- readRDS("rds/scRNA_annoted.rds")

library(Seurat)
Idents(scRNA)
table(scRNA$celltype)

# 如果idents是Seurat clusters，先将Idents修改为celltype
Idents(scRNA) <- "celltype"

# 提取小胶质细胞亚群
microglia_cells <- subset(scRNA, idents = "Microglia")

# 标准预处理流程
microglia_cells <- NormalizeData(microglia_cells)
microglia_cells <- FindVariableFeatures(microglia_cells)
microglia_cells <- ScaleData(microglia_cells)
microglia_cells <- RunPCA(microglia_cells)

# 使用UMAP降维并进行聚类，调整resolution
microglia_cells <- RunUMAP(microglia_cells, dims = 1:20)
microglia_cells <- FindNeighbors(microglia_cells, dims = 1:20)
microglia_cells <- FindClusters(microglia_cells, resolution = 0.2)  # 分辨率可以根据需求调整

table(microglia_cells$seurat_clusters)
# 使用ROGUE评估分辨率
library(Seurat)
library(tidyverse)
library(ROGUE)

# 提取表达矩阵
expr <- LayerData(microglia_cells, layer = "counts") %>% as.matrix()

# 过滤矩阵，去除在少于10个细胞中表达的基因
expr <- matr.filter(expr, min.cells = 10)

# 计算样本的 Shannon Entropy
ent.res <- SE_fun(expr)

# 绘制 Shannon Entropy 图
SEplot(ent.res)

# 计算 ROGUE 值
rogue.res <- rogue(expr, samples = "orig.ident",
                   labels = microglia_cells$seurat_clusters, platform = "UMI", span = 0.6)

# 查看 ROGUE 值
print(rogue.res)

# 绘制 ROGUE boxplot
rogue.boxplot(rogue.res)


# 可视化clustree
library(clustree)
clustree(microglia_cells)

# UMAP图展示小胶质细胞亚群
DimPlot(microglia_cells, reduction = "umap", label = TRUE) 
  + ggtitle("Microglia Subclusters")

library(ggrastr); library(Nebulosa)
ggrastr::rasterize(Nebulosa::plot_density(microglia_cells,
                       c("Mertk","Pparg"),
                       size = 0.2), dpi = 300)
FeaturePlot(microglia_cells, c("Mertk","Pparg"))

# 按样本 `orig.ident` 分组的 UMAP 图
DimPlot(microglia_cells, reduction = "umap", ncol = 3,
        split.by = "orig.ident") + ggtitle("UMAP of Microglia Subclusters by Sample")


# 找出每个小胶质细胞亚群的特异性marker
microglia_markers <- FindAllMarkers(microglia_cells, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# 保存marker到CSV文件
write.csv(microglia_markers, file = "microglia_markers.csv")
