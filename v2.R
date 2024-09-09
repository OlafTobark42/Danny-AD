library(Seurat)
library(limma)
library(patchwork) 
library(tidyverse)
library(SingleR)
BiocManager::install("singleR")
library(devtools)
devtools::install_github("dviraran/SingleR")
devtools::install_local("~/Downloads/SingleR-master-2.zip")
library("celldex")
# install.packages("BiocManager")
BiocManager::install("celldex")
library(SingleR)
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
  scRNAlist[[i]] <- RenameCells(scRNAlist[[i]], add.cell.id = names(ds.dir)[i])   
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
# Save file
if (dir.exists("rds")) saveRDS(scRNA,"rds/scRNA-raw-merge.rds") else {
  dir.create("rds")
  saveRDS(scRNA,"rds/scRNA-raw-merge.rds")
}

# Clean Up scRNAlist
rm(scRNAlist)
gc()

head(scRNA)
table(scRNA$orig.ident)

##==数据质控==#
# scRNA <- scRNA  #以后的分析使用整合的数据进行
##meta.data添加信息
proj_name <- data.frame(proj_name=rep("demo2",ncol(scRNA)))
rownames(proj_name) <- row.names(scRNA@meta.data)
scRNA <- AddMetaData(scRNA, proj_name)

##切换数据集
DefaultAssay(scRNA) <- "RNA"

##计算线粒体和红细胞基因比例
scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^MT-")
#计算红细胞比例
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB_m <- match(HB.genes, rownames(scRNA@assays$RNA)) 
HB.genes <- rownames(scRNA@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
scRNA[["percent.HB"]]<-PercentageFeatureSet(scRNA, features=HB.genes) 
#head(scRNA@meta.data)
col.num <- length(levels(as.factor(scRNA@meta.data$orig.ident)))

##绘制小提琴图
#所有样本一个小提琴图用group.by="proj_name"，每个样本一个小提琴图用group.by="orig.ident"
violin <-VlnPlot(scRNA, group.by = "proj_name",  
                 features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"), 
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
cl_tree <- clustree(scRNA)
ggsave("visual/clustree.pdf", width = 16, height = 9)

Idents(scRNA) <- scRNA$RNA_snn_res.0.4

scRNA <- RunUMAP(scRNA, reduction = "harmony", dims = 1:25)

table(scRNA$RNA_snn_res.0.4)
table(scRNA$seurat_clusters)
table(scRNA$orig.ident)
saveRDS(scRNA,"rds/scRNA-harmony-cluster.rds")

#选取主成分
# pc.num=1:25


##非线性降维
# scRNA <- RunTSNE(scRNA, dims = pc.num)
# scRNA <- RunUMAP(scRNA, dims = pc.num)

# DimPlot(scRNA, reduction = "tsne", label=T) 
# tsneplot <- DimPlot(scRNA, reduction = "tsne", group.by='orig.ident') 

DimPlot(scRNA, reduction = "umap", label=T) 
DimPlot(scRNA, reduction = "umap", group.by='orig.ident')

###definaion 细胞注释
##==鉴定细胞类型==##
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
library(BiocManager)
##之前安装过BiocManager，加载就可以。
BiocManager::install("SingleR", force = T)
library(SingleR)
BiocManager::install("celldex", force = T)
BiocManager::install("BiocSingular",force = T)
library(BiocSingular)
library(SingleR)
library("celldex")
#library(Seurat)

#dir.create("CellType")
#look cellex , first mouse cell 粗筛, then immune data label=label.fine
refdata <- celldex::ImmGenData()
# refdata <- MouseRNAseqData()
# testdata <- GetAssayData(scRNA, slot="data")
scRNA <- JoinLayers(scRNA)
testdata <- LayerData(scRNA)
clusters <- scRNA@meta.data$RNA_snn_res.0.4
saveRDS(scRNA, "rds/joinlayerd.rds")
scRNA$seurat_clusters <- scRNA$RNA_snn_res.0.4
#clusters
library(SingleR)
pred.mus <- SingleR(test = testdata, ref = refdata, labels = refdata$label.main, 
                    clusters = clusters, assay.type.test = "logcounts", assay.type.ref = "logcounts") #参数里没有test = 而是sc_data = 但升级2.2之后还是用test=
#pred.mus
celltype = data.frame(ClusterID=rownames(pred.mus), celltype=pred.mus$labels, stringsAsFactors = F)
dir.create("./CellType")
write.csv(celltype,"CellType/celltype.csv",row.names = F)
scRNA@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  scRNA@meta.data[which(scRNA@meta.data$seurat_clusters == celltype$ClusterID[i]),
                   'celltype'] <- celltype$celltype[i]}
p1 = DimPlot(scRNA, group.by="celltype", repel=T, label=T, label.size=5, reduction='tsne')
p2 = DimPlot(scRNA, group.by="celltype", repel=T, label=T, label.size=5, reduction='umap')
#p3 = p1+p2+ plot_layout(guides = 'collect')
dir.create("./visual")
library(ggplot2)
ggsave("./visual/umap_ImmGen.pdf",p2)


p <- ggplot(utest,aes(x= UMAP_1 , y = UMAP_2 ,color = cell_type)) +  
  geom_point(size = 1 , alpha =1 )  +  scale_color_manual(values = allcolour)


# 手动注释 microglia
scRNA@meta.data$celltype[which(scRNA@meta.data$seurat_clusters %in% c(2, 3, 7, 17, 20))] <- "Microglia"

# 手动注释 astrocytes
scRNA@meta.data$celltype[which(scRNA@meta.data$seurat_clusters %in% c(5, 23, 27))] <- "astrocytes"

# 手动注释 oligodendrocytes
scRNA@meta.data$celltype[which(scRNA@meta.data$seurat_clusters %in% c(0, 6, 12, 9))] <- "oligodendrocytes"
##保存数据
saveRDS(scRNA,'scRNA_annoted.rds')
scRNA <- readRDS('scRNA_annoted.rds')

FeaturePlot(scRNA, c("Cx3cr1","Ntsr2","Cldn11","Slc17a7","Cd3d","Cd3e","S100a9","Lyz2","Cd68","P2ry12","Acta2","Col1a1","Cldn5"))

## Cell Proportion
cellprop <- table(scRNA$orig.ident,scRNA$celltype)
cellprop <- proportions(cellprop,1)
as.data.frame(cellprop)
cellprop
proportions(cellprop,1)
colourCount = length(unique(cellprop$Var1))
ggplot(as.data.frame(cellprop))+
  geom_bar(aes(x =Var1, y= Freq, fill = Var2),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+ 
  labs(x='Sample',y = 'Ratio')



