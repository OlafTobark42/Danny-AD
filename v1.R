library(Seurat)
library(limma)
library(patchwork) 
library(tidyverse)
library(SingleR)
library("celldex")
library(ggplot2)

getwd()

fs=list.files('./rawdata',full.names=T)


scRNAlist = lapply(fs,function(folder){ 
  CreateSeuratObject(counts = Read10X(folder), 
                     project = folder )
})

#str(scRNAlist)
dim(scRNAlist[[2]]) #基因数  细胞数
table(scRNAlist[[2]]@meta.data$orig.ident)

##merge
scRNA1 <- merge(scRNAlist[[1]],scRNAlist[2:length(scRNAlist)])
table(scRNA1$orig.ident)
dim(scRNA1)
saveRDS(scRNA1,"mergedraws.rds")

scRNA1 <- readRDS("mergedraws.rds")

##==数据质控==#
scRNA <- scRNA1  #以后的分析使用整合的数据进行
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

scRNA1 <- scRNA

scRNA1 <- NormalizeData(scRNA1)
scRNA1 <- FindVariableFeatures(scRNA1, selection.method = "vst")
scRNA1 <- ScaleData(scRNA1, features = VariableFeatures(scRNA1))
scRNA1 <- RunPCA(scRNA1, features = VariableFeatures(scRNA1))
plot1 <- DimPlot(scRNA1, reduction = "pca", group.by="orig.ident")
plot2 <- ElbowPlot(scRNA1, ndims=30, reduction="pca") 
plotc <- plot1+plot2

#选取主成分
pc.num=1:25

##细胞聚类
scRNA1 <- FindNeighbors(scRNA1, dims = pc.num) 
scRNA1 <- FindClusters(scRNA1, resolution = 0.5)
table(scRNA1@meta.data$seurat_clusters)

##非线性降维
scRNA1 <- RunTSNE(scRNA1, dims = pc.num)
scRNA1 <- RunUMAP(scRNA1, dims = pc.num)

DimPlot(scRNA1, reduction = "tsne", label=T) 
tsneplot <- DimPlot(scRNA1, reduction = "tsne", group.by='orig.ident') 

DimPlot(scRNA1, reduction = "umap", label=T) 
umapplot <- DimPlot(scRNA1, reduction = "umap", group.by='orig.ident')

#合并tSNE与UMAP
plotc <- tsneplot+umapplot+ plot_layout(guides = 'collect')
pdf("umap+tsne1.pdf")
plotc
dev.off()

saveRDS(scRNA1, "harmonyed.rds")
scRNA1 <- readRDS("harmonyed.rds")
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
refdata <- ImmGenData()
testdata <- GetAssayData(scRNA1, slot="data")
clusters <- scRNA1@meta.data$seurat_clusters
#clusters

pred.mus <- SingleR(test = testdata, ref = refdata, labels = refdata$label.main, 
                    clusters = clusters, assay.type.test = "logcounts", assay.type.ref = "logcounts")
#pred.mus
celltype = data.frame(ClusterID=rownames(pred.mus), celltype=pred.mus$labels, stringsAsFactors = F)
dir.create("./CellType")
write.csv(celltype,"CellType/celltype.csv",row.names = F)
scRNA1@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  scRNA1@meta.data[which(scRNA1@meta.data$seurat_clusters == celltype$ClusterID[i]),
                   'celltype'] <- celltype$celltype[i]}
p1 = DimPlot(scRNA1, group.by="celltype", repel=T, label=T, label.size=5, reduction='tsne')
p2 = DimPlot(scRNA1, group.by="celltype", repel=T, label=T, label.size=5, reduction='umap')
#p3 = p1+p2+ plot_layout(guides = 'collect')
dir.create("./visual")
library(ggplot2)
ggsave("./visual/umap_ImmGen.pdf",p2)


p <- ggplot(utest,aes(x= UMAP_1 , y = UMAP_2 ,color = cell_type)) +  
  geom_point(size = 1 , alpha =1 )  +  scale_color_manual(values = allcolour)

##保存数据
saveRDS(scRNA1,'scRNA_annoted.rds')
scRNA1 <- readRDS('scRNA_annoted.rds')



## Cell Proportion
cellprop <- table(scRNA1$orig.ident,scRNA1$celltype)
cellprop <- proportions(cellprop,1)
as.data.frame(cellprop)
cellprop
proportions(cellprop,1)
colourCount = length(unique(cellprop$Var1))
ggplot(as.data.frame(cellprop))+
  geom_bar(aes(x =Var1, y= Freq, fill = Var2),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+ 
  labs(x='Sample',y = 'Ratio')



##find all markers
dir.create("cell_identify")
diff.wilcox = FindAllMarkers(scRNA1)
all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val<0.05)
top50 = all.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
write.csv(top50, "cell_identify/top50_diff_genes_wilcox.csv", row.names = F)
write.csv(all.markers,"cell_identify/all_markers_diff_genes_wilcox.csv",row.names=F)


#cut ""
vector(top50[top50$cluster==3,1])
#from error a <- '"gene1","gene2","gene3" ';b=gsub('["]', '', a)
#b to CellMarker 2.0 http://bio-bigdata.hrbmu.edu.cn/CellMarker/CellMarker_annotation.jsp
a <- '"Sesn1", "Tsc22d3", "Zfp36l2", "Cirbp", "Apobec1", "Aldh2", "Rbm3", "Serinc3", "Apoe", 
"Snx5", "Btla", "Ddit4", "H2-Aa", "Birc3", "Ppp1r21", "Pik3ip1", "Neurl3", "Pkig", "Cd79b", 
"Fam107b", "Cyp4f18", "Dok3", "Ptpn22", "H2-Eb1", "H2-Oa", "Fkbp5", "Serpinb1a", "Sypl", 
"Cd200", "Cd79a", "H2-Ob", "Vars", "Ston1", "Scd1", "Pou2af1", "H3f3b", "Abca1", "Ebf1", 
"Ighd", "Hvcn1", "H2-Ab1", "Cd74", "Gm8369", "Mylip", "1110059E24Rik", "Sgk3", "Btk", 
"Haao", "Trp53i11", "Herc4" ';b=gsub('["]', '', a)
b

###
markers_df <- FindMarkers(object = scRNA1, ident.1 = 9, min.pct = 0.25)
markers_genes =  head(x = top50_diff_genes_wilcox$gene, n = 5)
print(x = head(markers_df))
markers_genes =  rownames(head(x = markers_df, n = 5))
VlnPlot(object = scRNA1, features =markers_genes,log =T )
FeaturePlot(object = scRNA1, features=markers_genes )
FeaturePlot(scRNA1,features = markers_genes)

# 官方示例
cell_typeA_marker_gene_list <- list(c("Gene1", "Gene2", "Gene3", "Gene4"))
object <- AddModuleScore(object = object, features = cell_typeA_marker_gene_list, name = "cell_typeA_score")
FeaturePlot(object = object, features = "cell_typeA_score1")

# 给出一段GitHub文献中的代码，学习下大佬的使用方法
b_cell <- c("MS4A1") # ref5
macrophage <- c("CD68", "IDO1") # ref5
plasmacytoid_dendritic_cell <- c("CLEC4C","NRP1") # ref5
erythrocyte <- c("HBB", "HBA1", "HBA2") # ref6
cytotoxic_t_cell <- c("GZMA", "GZMK", "IFNG") # ref5
regulatory_t_cell <- c("FOXP3", "IL2RA", "IKZF2") # ref5
t_helper <- c("CXCL13","PDCD1","FABP5") # ref5
naive_t_cell <- c("CCR7", "IL7R", "LEF1") # ref5
progenitor <- c("CD34") # ref5
mast_cell <- c("TPSAB1","TPSB2", "KIT", "GATA1", "GATA2") # ref7

# marker基因打分函数
score_marker_genes <- function(seurat_object, nbins=24){
  AddModuleScore(seurat_object, 
                 features = list(b_cell, macrophage, plasmacytoid_dendritic_cell, erythrocyte, cytotoxic_t_cell, regulatory_t_cell, t_helper, naive_t_cell, progenitor, mast_cell),
                 name=c("b_cell", "macrophage", "plasmacytoid_dendritic_cell", "erythrocyte", "cytotoxic_t_cell", "regulatory_t_cell", "t_helper", "naive_t_cell", "progenitor", "mast_cell"),
                 nbin=nbins)
}
sc <- score_marker_genes(scRNA1)

# marker基因可视化函数
plot_marker_genes <- function(seurat_object, image_name){
  ggsave(file = str_glue('../figures/{image_name}.pdf'),
         plot = FeaturePlot(seurat_object,
                            features = c("b_cell1", "macrophage2", "plasmacytoid_dendritic_cell3", "erythrocyte4", "cytotoxic_t_cell5", "regulatory_t_cell6", "t_helper7", "naive_t_cell8", "progenitor9", "mast_cell10"),
                            min.cutoff = "q10", max.cutoff = "q90",
                            ncol=4, label=TRUE, order = TRUE),
         device = "pdf", width = 50, height = 40, units = "cm")
}
plot_marker_genes(hl, "hl_markers")

cgs = list(
  Epi = c("EPCAM","PAX8","KRT18","CD24","KRT19","SCGB2A2","KRT5","KRT15"),
  Meyloid = c("CD68","LYZ","MARCO","AIF1","TYROBR","MS4A6A","CD1E","IL3RA","LAMP3"),
  T_cell = c("CD3D",'CD3E','TRAG','CD3G','CD2'),
  B_cell = c("CD79A","CD79B","IGKC","CD19","MZB1","MS4A1"),
  Endo = c("CLDN5","PECAM1","VWF","FLT1","RAMP2"),
  Fibro = c("COL1A1","COL1A2","COL3A1","BGN","DCN","POSTN","C1R")
)
testlist <- list(regulatory_t_cell =c("FOXP3", "IL2RA", "IKZF2"), # ref5
                 t_helper =c("CXCL13","PDCD1","FABP5") ,# ref5
                 naive_t_cell =c("CCR7", "IL7R", "LEF1"), # ref5
                 progenitor =c("CD34"), # ref5
                 mast_cell =c("TPSAB1","TPSB2", "KIT", "GATA1", "GATA2") # ref7
)

p_umap <- DimPlot(scRNA1, reduction = "umap", group.by = "celltype", label = T)
p <- DotPlot(scRNA1, features = testlist) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + p_umap + plot_layout(widths = c(2, 1))
save_plot("Dotplot_Celltype_UMAP_clustering.png", p, base_height = 8, base_aspect_ratio = 2.4, base_width = NULL, dpi=600)
save_plot("Dotplot_Celltype_UMAP_clustering.pdf", p, base_height = 8, base_aspect_ratio = 2.4, base_width = NULL)



