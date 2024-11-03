rm(list = ls()); gc()
library(Seurat)
library(limma)
library(patchwork) 
library(tidyverse)
library(SingleR)
# BiocManager::install("singleR")
library(devtools)
# devtools::install_github("dviraran/SingleR")
# devtools::install_local("~/Downloads/SingleR-master-2.zip")
library("celldex")
# BiocManager::install("celldex")
library(ggplot2)

getwd()

ds.dir=list.files('./rawdata',full.names=T)
names(ds.dir) <- list.files('./rawdata')
# First Specify the clasee of scRNAlist
scRNAlist <- list()
for (i in 1:length(ds.dir)) {
  scRNAlist[[i]] <- CreateSeuratObject(counts=Read10X(data.dir = ds.dir[i]), 
                                       project=names(ds.dir)[i], min.cells=3, min.features = 200)
  # 给细胞barcode加个前缀，防止合并后barcode重名
  # scRNAlist[[i]] <- RenameCells(scRNAlist[[i]], add.cell.id = names(ds.dir)[i])   
  # 计算线粒体基因比例
  if (T) {    
    scRNAlist[[i]][["percent.mt"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern = "^mt-")   # Hm: MT
  }
  # 计算核糖体基因比例
  if (T) {
    scRNAlist[[i]][["percent.rb"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern = "^Rp[sl]")  # Hm: RP[SL]
  }
  #计算红细胞基因比例
  if (T) {
    library(stringr)
    HB.genes <- str_to_title(tolower(c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")))
    HB.genes <- CaseMatch(HB.genes, rownames(scRNAlist[[i]]))
    scRNAlist[[i]][["percent.HB"]]<-PercentageFeatureSet(scRNAlist[[i]], features=HB.genes) 
  }
  # this subset step should be in 02-qc
  # scRNAlist[[i]] <- subset(scRNAlist[[i]], subset = percent.mt < 10) 
}

# Merge multi-samples into one object
scRNA <- merge(scRNAlist[[1]],scRNAlist[2:length(scRNAlist)])
table(scRNA$orig.ident)
head(scRNA)
# Clean Up scRNAlist
rm(scRNAlist); gc()

##==数据质控==#
##meta.data添加信息
proj_name <- data.frame(proj_name=rep("danny-AD",ncol(scRNA)))
rownames(proj_name) <- row.names(scRNA@meta.data)
scRNA <- AddMetaData(scRNA, proj_name)

##切换数据集
DefaultAssay(scRNA) <- "RNA"

col.num <- length(levels(as.factor(scRNA@meta.data$orig.ident)))

##绘制小提琴图
#所有样本一个小提琴图用group.by="proj_name"，每个样本一个小提琴图用group.by="orig.ident"
violin <-VlnPlot(scRNA, group.by = "proj_name",  
                 features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb","percent.HB"), 
                 cols =rainbow(col.num), 
                 pt.size = 0.01, #不需要显示点，可以设置pt.size = 0
                 ncol = 4) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) 
dir.create("./QC")
ggsave("QC/vlnplot_before_qc.pdf", plot = violin, width = 12, height = 6) 
ggsave("QC/vlnplot_before_qc.png", plot = violin, width = 12, height = 6)  
plot1 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "percent.HB")
pearplot <- CombinePlots(plots = list(plot1, plot2, plot3), nrow=1, legend="none") 
ggsave("QC/pearplot_before_qc.pdf", plot = pearplot, width = 12, height = 5) 
ggsave("QC/pearplot_before_qc.png", plot = pearplot, width = 12, height = 5)

##设置质控标准
print(c("请输入允许基因数和核糖体比例，示例如下：", "minGene=500", "maxGene=4000", "pctMT=20"))
minGene=500
maxGene=3000
pctMT=10

##数据质控
scRNA <- subset(scRNA, subset = nFeature_RNA > minGene & nFeature_RNA < maxGene & percent.mt < pctMT)
col.num <- length(levels(as.factor(scRNA@meta.data$orig.ident)))
violin <-VlnPlot(scRNA, group.by = "proj_name",
                 features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"), 
                 cols =rainbow(col.num), 
                 pt.size = 0.1, 
                 ncol = 4) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) 
ggsave("QC/vlnplot_after_qc.pdf", plot = violin, width = 12, height = 6) 
ggsave("QC/vlnplot_after_qc.png", plot = violin, width = 12, height = 6)


dim(scRNA)   #查看基因数和细胞总数
table(scRNA@meta.data$orig.ident)  #查看每个样本的细胞数

library(harmony)
scRNA <- NormalizeData(scRNA) %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA(verbose=FALSE)
scRNA <- RunHarmony(scRNA, group.by.vars = "orig.ident")
scRNA <- FindNeighbors(scRNA, reduction = "harmony", dims = 1:25) %>% 
  FindClusters(resolution = seq(0.1, 1, by = 0.1))

library(clustree)
# cl_tree <- clustree(scRNA)
ggsave("visual/clustree.pdf", width = 16, height = 9)

Idents(scRNA) <- scRNA$RNA_snn_res.0.5
scRNA$seurat_clusters <- scRNA$RNA_snn_res.0.5

library(dplyr)
# 提取 meta.data 并去掉 RNA_snn_res.0.1-1 的列
scRNA@meta.data <- scRNA@meta.data %>%
  select(-starts_with("RNA_snn_res."))
head(scRNA)

# 创建分组变量并按因子排序
scRNA$group <- ifelse(scRNA$orig.ident %in% c("a1", "a2", "a3"), "AD + WT",
                                ifelse(scRNA$orig.ident %in% c("c1", "c2", "c3"), "Control",
                                       "AD + Mertk-/-"))

# 将分组按指定顺序设为因子
scRNA$group <- factor(scRNA$group, levels = c("Control", "AD + WT", "AD + Mertk-/-"))

# 检查新的分组变量
table(scRNA$group)

# UMAP降维
scRNA <- RunUMAP(scRNA, reduction = "harmony", dims = 1:25)

DimPlot(scRNA, reduction = "umap", label=T) 
DimPlot(scRNA, reduction = "umap", split.by = 'group')

##==鉴定细胞类型==##
FeaturePlot(scRNA, brain_gene[["Oligodendrocyte progenitor cells"]])
FeaturePlot(scRNA, immune_gene[["Neutrophils"]])
FeaturePlot(scRNA, brain_gene[["Vascular cells"]])
FeaturePlot(scRNA, brain_gene[["Endothelial vascular cells"]])

scRNA <- JoinLayers(scRNA)
library(devtools)
devtools::install_local("~/Downloads/presto-master.zip")
marker_11 <- FindMarkers(scRNA, ident.1 = 11)
head(marker_11, 20)

scRNA@meta.data$celltype <- NA
scRNA@meta.data$celltype[which(scRNA@meta.data$seurat_clusters %in% c(25))] <- "T cells"
scRNA@meta.data$celltype[which(scRNA@meta.data$seurat_clusters %in% c(1,6,10))] <- "Endothelial Cells"
scRNA@meta.data$celltype[which(scRNA@meta.data$seurat_clusters %in% c(13))] <- "VSMCs"
scRNA@meta.data$celltype[which(scRNA@meta.data$seurat_clusters %in% c(11))] <- "Ependymal Cells"

# neuron
scRNA@meta.data$celltype[which(scRNA@meta.data$seurat_clusters %in% c(8,23,24,27,28))] <- "Neuron"

# 手动注释 microglia
scRNA@meta.data$celltype[which(scRNA@meta.data$seurat_clusters %in% c(2,5,4,12,17))] <- "Microglia"
scRNA@meta.data$celltype[which(scRNA@meta.data$seurat_clusters %in% c(9,19,29))] <- "Macrophages"

# 手动注释 astrocytes
scRNA@meta.data$celltype[which(scRNA@meta.data$seurat_clusters %in% c(7,14,18))] <- "Astrocytes"

# 手动注释 oligodendrocytes
scRNA@meta.data$celltype[which(scRNA@meta.data$seurat_clusters %in% c(0,3,15,16,26))] <- "Oligodendrocytes"

# show results
table(scRNA$celltype)
##保存数据
saveRDS(scRNA,'scRNA_annoted.rds')
scRNA <- readRDS('scRNA_annoted.rds')

FeaturePlot(scRNA, c("Cx3cr1","Ntsr2","Cldn11","Slc17a7","Cd3d","Cd3e","S100a9","Lyz2","Cd68","P2ry12","Acta2","Col1a1","Cldn5"))

## Cell Proportion
cellprop <- table(scRNA$group,scRNA$celltype)
cellprop <- proportions(cellprop,1)
cellprop <- as.data.frame(cellprop)

colourCount = length(unique(cellprop$Var1))
ggplot(as.data.frame(cellprop))+
  geom_bar(aes(x =Var1, y= Freq, fill = Var2),stat = "identity",width = 0.7,linewidth = 0.5,colour = '#222222')+ 
  labs(x='Sample',y = 'Ratio')

saveRDS(scRNA, "rds/manual_annoted.rds")

