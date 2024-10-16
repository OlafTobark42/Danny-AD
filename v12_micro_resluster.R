library(Seurat)

seu <- readRDS("~/Downloads/scAD/rds/done_annotated.rds")
head(seu)
DimPlot(seu, label = T, group.by = "celltype")

FeaturePlot(seu, "Pros1", split.by = "group")
VlnPlot(seu, "Gas6", group.by = "celltype", split.by = "group")



micro_sub <- readRDS("~/Downloads/scAD/rds/micro_sub.rds")
DimPlot(micro_sub, group.by = "seurat_clusters", label = T)
FeaturePlot(micro_sub, "Gas6")
VlnPlot(micro_sub, "Gas6", group.by = "seurat_clusters", split.by = "group")

VlnPlot(micro_sub, "Pros1", pt.size = 0,
        group.by = "seurat_clusters", split.by = "group")

library(tidyverse)
table(seu$celltype)
neuron_sub <- subset(seu, celltype %in% "Neuron")
VlnPlot(neuron_sub, "Pros1", pt.size = 0.5, group.by = "group")



micro4 <- subset(seu, celltype %in% "Microglia")
  # 标准预处理流程
micro4 <- NormalizeData(micro4)
micro4 <- FindVariableFeatures(micro4)
micro4 <- ScaleData(micro4)
micro4 <- RunPCA(micro4)

# 使用UMAP降维并进行聚类，调整resolution
micro4 <- RunUMAP(micro4, dims = 1:20)
micro4 <- FindNeighbors(micro4, dims = 1:20)
micro4 <- FindClusters(micro4, resolution = 1)  # 分辨率可以根据需求调整

table(micro4$seurat_clusters)
DimPlot(micro4, label = T)
VlnPlot(micro4, "Mertk", pt.size = 0,
        group.by = "seurat_clusters", split.by = "group")

FeaturePlot(micro4, "Cd68", split.by = "group")
VlnPlot(micro4, "Cd68", 
        group.by = "seurat_clusters", split.by = "group")

Idents(micro_sub) <- micro_sub$group
table(micro_sub$seurat_clusters)
m3 <- subset(micro_sub, seurat_clusters %in% 3)
deg3 <- FindMarkers(m3, ident.1 = "AD + WT", ident.2 = "Control")
deg3["Pparg",]
