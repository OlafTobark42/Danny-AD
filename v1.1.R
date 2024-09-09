scRNA1 <- readRDS("harmonyed.rds")

scRNA1$orig.ident[scRNA1$orig.ident == "a1",]
# Assuming scRNA1$orig.ident is a dataframe or a column from a dataframe
scRNA1$group <- NA
scRNA1$group <- ifelse(scRNA1$orig.ident %in% c('./rawdata/c1', 
                                                './rawdata/c2', 
                                                './rawdata/c3'), 'cont', scRNA1$group)

table(scRNA1$group)
library(Seurat)
DimPlot(scRNA1, split.by = "group")
FeaturePlot(scRNA1, "Pparg", split.by = "group", max.cutoff = 1)

ko <- subset(scRNA1, group == "ko")
cont <- subset(scRNA1, group == "cont")
mertk <- subset(scRNA1, seurat_clusters == "2")
Idents(scRNA1) <- scRNA1$seurat_clusters
scRNA1 <- JoinLayers(scRNA1)
FindMarkers(scRNA1, ident.1 = "2")
head(scRNA1)

FeaturePlot(scRNA1, "Cx3cr1")
FeaturePlot(scRNA1, "P2ry12")
FeaturePlot(scRNA1, "Tmem119")
FeaturePlot(scRNA1, "Cd68")
FeaturePlot(scRNA1, "Ms4a1")
