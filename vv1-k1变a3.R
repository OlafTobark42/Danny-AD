rm(list = ls()); gc()
library(Seurat)
library(limma)
library(patchwork) 
library(tidyverse)
library(SingleR)
# BiocManager::install("SingleR")
library(devtools)
# devtools::install_github("dviraran/SingleR")
# devtools::install_local("~/Downloads/SingleR-master-2.zip")
library("celldex")
# BiocManager::install("celldex")
library(ggplot2);library(ggsci); library(RColorBrewer);library(ggalluvial)

getwd()

ds.dir=list.files('./rawdata',full.names=T)
names(ds.dir) <- list.files('./rawdata')
# First Specify the clasee of seulist
seulist <- list()
for (i in 1:length(ds.dir)) {
  seulist[[i]] <- CreateSeuratObject(counts=Read10X(data.dir = ds.dir[i]), 
                                       project=names(ds.dir)[i], min.cells=3, min.features = 200)
  # 给细胞barcode加个前缀，防止合并后barcode重名
  # seulist[[i]] <- RenameCells(seulist[[i]], add.cell.id = names(ds.dir)[i])   
  # 计算线粒体基因比例
  if (T) {    
    seulist[[i]][["Percent.Mito"]] <- PercentageFeatureSet(seulist[[i]], pattern = "^mt-")   # Hm: MT
  }
  # 计算核糖体基因比例
  if (T) {
    seulist[[i]][["Percent.Ribo"]] <- PercentageFeatureSet(seulist[[i]], pattern = "^Rp[sl]")  # Hm: RP[SL]
  }
  #计算红细胞基因比例
  if (F) {
    library(stringr)
    HB.genes <- str_to_title(tolower(c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")))
    # HB.genes <- CaseMatch(HB.genes, rownames(seulist[[i]]))
    seulist[[i]][["percent.HB"]]<-PercentageFeatureSet(seulist[[i]], features=HB.genes) 
  }
  # this subset step should be in 02-qc
  # seulist[[i]] <- subset(seulist[[i]], subset = percent.mt < 10) 
}

# Merge multi-samples into one object
seu <- merge(seulist[[1]],seulist[2:length(seulist)])
table(seu$orig.ident)
head(seu)
# Clean Up seulist
rm(seulist); gc()

seu$orig.ident <- factor(seu$orig.ident, labels = (c("Ctrl1","Ctrl2","Ctrl3","AßO1","AßO2","AßO3")),
                           levels = c("c1","c2","c3","a1","a2","a3"))

##==数据质控==#
# ##meta.data添加信息
# proj_name <- data.frame(proj_name=rep("danny-AD",ncol(seu)))
# rownames(proj_name) <- row.names(seu@meta.data)
# seu <- AddMetaData(seu, proj_name)

##切换数据集
DefaultAssay(seu) <- "RNA"

col.num <- length(levels(as.factor(seu@meta.data$orig.ident)))
library(ggsci);library(RColorBrewer)
##绘制小提琴图
#所有样本一个小提琴图用group.by="proj_name"，每个样本一个小提琴图用group.by="orig.ident"
violin <-VlnPlot(seu, group.by = "orig.ident",  
                 features = c("nFeature_RNA", "nCount_RNA", "Percent.Mito", "Percent.Ribo"), 
                 # cols = ggsci::pal_npg(palette = "nrc")(col.num), 
                 col = RColorBrewer::brewer.pal(col.num, "Set2"),
                 pt.size = 0, #不需要显示点，可以设置pt.size = 0
                 ncol = 4) + 
  scale_y_continuous(limits = c(0.01, NA)) +  # 设置 Y 轴起点略高于 0
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) 
violin

dir.create("./QC")
ggsave("QC/vlnplot_before_qc.pdf", plot = violin, width = 12, height = 6) 

dim(seu)   #查看基因数和细胞总数
table(seu@meta.data$orig.ident)  #查看每个样本的细胞数

##设置质控标准
print(c("请输入允许基因数和核糖体比例，示例如下：", "minGene=500", "maxGene=4000", "pctMT=20"))
minGene=500
maxGene=3000
pctMT=10

##数据质控
seu <- subset(seu, subset = nFeature_RNA > minGene & nFeature_RNA < maxGene & Percent.Mito < pctMT)
# col.num <- length(levels(as.factor(seu@meta.data$orig.ident)))
violin <-VlnPlot(seu, group.by = "orig.ident",  
                 features = c("nFeature_RNA", "nCount_RNA", "Percent.Mito", "Percent.Ribo"), 
                 # cols = ggsci::pal_npg(palette = "nrc")(col.num), 
                 col = RColorBrewer::brewer.pal(col.num, "Set2"),
                 pt.size = 0, #不需要显示点，可以设置pt.size = 0
                 ncol = 4) + 
  scale_y_continuous(limits = c(0.01, NA)) +  # 设置 Y 轴起点略高于 0
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) 
violin
ggsave("QC/vlnplot_after_qc.pdf", plot = violin, width = 12, height = 6) 
# ggsave("QC/vlnplot_after_qc.png", plot = violin, width = 12, height = 6)


dim(seu)   #查看基因数和细胞总数
table(seu@meta.data$orig.ident)  #查看每个样本的细胞数

library(harmony)
seu <- NormalizeData(seu) %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA(verbose=FALSE)
seu <- RunHarmony(seu, group.by.vars = "orig.ident")
seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:25) %>% 
  FindClusters(resolution = seq(0.1, 1, by = 0.1))

library(clustree)
library(ggplot2)

# 创建一个包含 30 种颜色的自定义调色板
custom_colors <- c(
  RColorBrewer::brewer.pal(12, "Paired"), 
  RColorBrewer::brewer.pal(8, "Set3"), 
  ggsci::pal_npg("nrc", alpha = 1)(10)
)

# 用于 ggplot 或 Seurat 图
scale_color_manual(values = custom_colors)

cl_tree <- clustree(seu) + RColorBrewer::brewer.pal(n = 10,name = "Paired")
cl_tree
ggsave("visual/clustree.pdf", cl_tree, width = 12, height = 12)

seu$seurat_clusters <- seu$RNA_snn_res.0.1
Idents(seu) <- seu$seurat_clusters

library(dplyr)
# 提取 meta.data 并去掉 RNA_snn_res.0.1-1 的列
seu@meta.data <- seu@meta.data %>%
  select(-starts_with(c("RNA_snn_res.","Percent.")))
head(seu)

# 创建分组变量并按因子排序
seu$group <- ifelse(seu$orig.ident %in% c("AßO1", "AßO2", "AßO3"), "AßO",
                      "Ctrl")

# 将分组按指定顺序设为因子
seu$group <- factor(seu$group, levels = c("Ctrl", "AßO"))

# 检查新的分组变量
table(seu$group)

# UMAP降维
seu <- RunUMAP(seu, reduction = "harmony", dims = 1:25)


# 创建一个包含 30 种颜色的自定义调色板
custom_colors <- c(
  RColorBrewer::brewer.pal(8, "Set1"), 
  RColorBrewer::brewer.pal(8, "Set2"), 
  ggsci::pal_npg("nrc", alpha = 1)(10)
)

# 用于 ggplot 或 Seurat 图
scale_color_manual(values = custom_colors)

DimPlot(seu, reduction = "umap", label=T) +scale_color_manual(values = custom_colors) 
ggsave("visual/umap-cluster.pdf", width = 8, height = 6)
DimPlot(seu, reduction = "umap", split.by = 'group')


# Doublet Find ------------------------------------------------------------
if (F) {
  library(Seurat)
  library(DoubletFinder)
  # BiocManager::install("scDblFinder")  似乎适用于sce格式
  # 3. 运行 DoubletFinder
  seu <- JoinLayers(seu)
  sweep.res.list <- paramSweep(seu, PCs = 1:10, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  # 提取出全局最优的pK值，储存于"pK_bcmvn"
  pK_bcmvn <- as.numeric(bcmvn$pK[which.max(bcmvn$BCmetric)]) 
  
  # 根据数据集确定期望双胞体数量
  annotations <- seu@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)
  nExp_poi <- round(0.075 * nrow(seu@meta.data))  # 假设 7.5% 的双胞体比例
  nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))  # 计算异源双细胞数量
  
  save.image(file = "cal_dbl.rda")
  # Error: vector memory limit of 32.0 Gb reached, see mem.maxVSize()
  seu <- doubletFinder(seu, PCs = 1:10, pN = 0.25, pK = pK_bcmvn, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
  
  # 4. 可视化双胞体分布
  DimPlot(seu, reduction = "umap", group.by = "DF.classifications_0.25_0.09_XX")
  
  # 5. 过滤双胞体
  seu <- subset(seu, subset = DF.classifications_0.25_0.09_XX == "Singlet")
}

# 换用另一个包 ------------------------------------------------------------------

library(Seurat)
library(SingleCellExperiment)
library(scDblFinder)

# 转换 Seurat 对象到 SingleCellExperiment 格式
sce <- as.SingleCellExperiment(seu)

# 指定样本来源
sce$sample <- seu$orig.ident  # 假设样本来源在 Seurat 的 orig.ident 中

# 运行 scDblFinder，指定样本
sce <- scDblFinder(sce, samples = sce$sample)

# 查看预测的双胞体分类
table(sce$scDblFinder.class)

# 将结果添加回 Seurat 对象
seu$doublet_class <- sce$scDblFinder.class

brewer.pal.info
RColorBrewer::brewer.pal(n = 3,"Greys")
# 可视化双胞体分布
dbl <- DimPlot(seu, group.by = "doublet_class",cols = RColorBrewer::brewer.pal(n = 3,"Greys"))
ggsave(filename = "visual/dbl.pdf", width = 8, height = 6)

# 保留非双胞体
seu <- subset(seu, subset = doublet_class == "singlet")


# Annotation --------------------------------------------------------------

source("code/markers.R")

##==鉴定细胞类型==##
common <- FeaturePlot(seu, c("Cx3cr1","Ntsr2","Cldn11","Slc17a7","Cd3d","Cd3e","S100a9","Ms4a1","Cd19",
                               "Lyz2","Cd68","P2ry12","Acta2","Col1a1","Cldn5"))
ggsave("visual/common_marker.pdf", common, width = 10 ,height = 14)
FeaturePlot(seu, brain_gene[["Oligodendrocyte progenitor cells"]])
FeaturePlot(seu, immune_gene[["Neutrophils"]])
FeaturePlot(seu, immune_gene[["DCs"]])
FeaturePlot(seu, brain_gene[["Vascular cells"]])
FeaturePlot(seu, brain_gene[["Neurons"]])
FeaturePlot(seu, brain_gene[["Oligodendrocytes"]])
FeaturePlot(seu, brain_gene[["Fibroblasts"]])
FeaturePlot(seu, brain_gene[["Choroid cells"]])
FeaturePlot(seu, marker_genes[["Monocytes"]])


seu <- JoinLayers(seu)
# library(devtools)
# devtools::install_local("~/Downloads/presto-master.zip")
marker_11 <- FindMarkers(seu, ident.1 = 11)
head(marker_11, 30)
marker_7 <- FindMarkers(seu, ident.1 = 7)
head(marker_7, 30)
marker_14 <- FindMarkers(seu, ident.1 = 14)
head(marker_14, 30)

seu@meta.data$celltype <- NA
seu@meta.data$celltype[which(seu@meta.data$seurat_clusters %in% c(12))] <- "Fibroblasts"
seu@meta.data$celltype[which(seu@meta.data$seurat_clusters %in% c(2,10))] <- "Endothelial Cells"
seu@meta.data$celltype[which(seu@meta.data$seurat_clusters %in% c(5))] <- "VSMCs"
seu@meta.data$celltype[which(seu@meta.data$seurat_clusters %in% c(11))] <- "Ependymal Cells"
seu@meta.data$celltype[which(seu@meta.data$seurat_clusters %in% c(7))] <- "Choroid Plexus Cells"

# neuron
seu@meta.data$celltype[which(seu@meta.data$seurat_clusters %in% c(4,14))] <- "Neuron"

# 手动注释 microglia
seu@meta.data$celltype[which(seu@meta.data$seurat_clusters %in% c(1,9))] <- "Microglia"
seu@meta.data$celltype[which(seu@meta.data$seurat_clusters %in% c(6))] <- "Macrophages"

# 手动注释 astrocytes
seu@meta.data$celltype[which(seu@meta.data$seurat_clusters %in% c(3,14))] <- "Astrocytes"

# 手动注释 oligodendrocytes
seu@meta.data$celltype[which(seu@meta.data$seurat_clusters %in% c(0,8))] <- "Oligodendrocytes"

# 手动 immune cells
seu@meta.data$celltype[which(seu@meta.data$seurat_clusters %in% c(13))] <- "Immune Cells"

# show results
table(seu$celltype)

seu$celltype <- factor(seu$celltype,
                         levels = c("Microglia", 
                                    setdiff(unique(seu$celltype),"Microglia")
                                    )
                         )
table(seu$celltype)
Idents(seu) <- seu$celltype
##保存数据
saveRDS(seu,'rds/seu_annoted.rds')
# seu <- readRDS('rds/seu_annoted.rds')

# umap-colorful -----------------------------------------------------------
# 加载必要的库
library(Seurat)
library(ggplot2)
library(ggsci)
library(RColorBrewer)

# 设置输出目录
output_dir <- "visual/umap-colorful"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# 获取所有细胞类型
celltypes <- unique(seu$celltype)
n_celltypes <- length(celltypes)

# 获取 ggsci 包中所有可用的配色方案
ggsci_palettes <- list(
  "ggsci_npg" = pal_npg("nrc"),
  "ggsci_aaas" = pal_aaas("default"),
  "ggsci_jco" = pal_jco("default"),
  "ggsci_lancet" = pal_lancet("lanonc"),
  "ggsci_ucscgb" = pal_ucscgb("default"),
  "ggsci_d3" = pal_d3("category20"),
  "ggsci_material_blue" = pal_material("blue"), # 示例颜色主题
  "ggsci_material_red" = pal_material("red"),   # 示例颜色主题
  "ggsci_futurama" = pal_futurama("planetexpress"),
  "ggsci_rickandmorty" = pal_rickandmorty("schwifty")
)
ggsci_palettes[[1]]

# 加载必要的库
library(RColorBrewer)

# 获取 RColorBrewer 中所有 "qual" 和 "div" 类型的配色方案
rcb_palettes <- rownames(brewer.pal.info[brewer.pal.info$category %in% c("qual", "div"), ])

# 获取每个配色方案的最大颜色数
rcb_palette_info <- brewer.pal.info[rownames(brewer.pal.info) %in% rcb_palettes, "maxcolors"]

# 检查结果
rcb_palettes
rcb_palette_info


# 创建补充颜色函数
get_filled_palette <- function(colors, n_needed) {
  # 如果颜色不足，补充随机颜色；否则返回前 n_needed 个颜色
  if (length(colors) < n_needed) {
    extra_colors <- grDevices::colors()[sample(1:length(grDevices::colors()), n_needed - length(colors))]
    c(colors, extra_colors)
  } else {
    colors[1:n_needed]
  }
}

# 遍历 ggsci 配色方案
for (palette_name in names(ggsci_palettes)) {
  palette_function <- ggsci_palettes[[palette_name]]
  
  # 获取当前配色方案颜色，补足颜色数量
  palette_colors <- get_filled_palette(palette_function(n_celltypes), n_celltypes)
  
  # 创建 DimPlot 并应用配色方案
  dimplot <- DimPlot(
    seu,
    group.by = "celltype",
    label = TRUE,       # 显示标签
    label.size = 4      # 标签字体大小
  ) +
    scale_color_manual(values = palette_colors) + # 应用配色方案
    theme_minimal() +  # 使用简单的背景主题
    theme(
      legend.position = "right",         # 图例显示在右侧
      panel.grid = element_blank(),      # 去掉网格线
      axis.text = element_blank(),       # 去掉坐标轴文字
      axis.ticks = element_blank(),      # 去掉坐标轴刻度
      axis.title = element_blank(),      # 去掉坐标轴标题
      plot.title = element_blank()       # 去掉 UMAP 图标题
    )
  
  # 保存为 PDF，文件名包含配色方案的名称
  output_file <- file.path(output_dir, paste0("DimPlot_", palette_name, ".pdf"))
  ggsave(
    filename = output_file,
    plot = dimplot,
    width = 8,  # 设置 PDF 宽度
    height = 6  # 设置 PDF 高度
  )
  
  message("已保存: ", output_file)
}

# 遍历 RColorBrewer 配色方案
for (palette_name in rcb_palettes) {
  max_colors <- brewer.pal.info[palette_name, "maxcolors"]
  
  # 如果配色方案颜色不足，则补充；否则使用完整配色
  if (max_colors < n_celltypes) {
    message("配色方案 ", palette_name, " 颜色不足，补充随机颜色。")
    palette_colors <- get_filled_palette(brewer.pal(max_colors, palette_name), n_celltypes)
  } else {
    palette_colors <- brewer.pal(n_celltypes, palette_name)
  }
  
  # 创建 DimPlot 并应用配色方案
  dimplot <- DimPlot(
    seu,
    group.by = "celltype",
    label = TRUE,       # 显示标签
    label.size = 4      # 标签字体大小
  ) +
    scale_color_manual(values = palette_colors) + # 应用配色方案
    theme_minimal() +  # 使用简单的背景主题
    theme(
      legend.position = "right",         # 图例显示在右侧
      panel.grid = element_blank(),      # 去掉网格线
      axis.text = element_blank(),       # 去掉坐标轴文字
      axis.ticks = element_blank(),      # 去掉坐标轴刻度
      axis.title = element_blank(),      # 去掉坐标轴标题
      plot.title = element_blank()       # 去掉 UMAP 图标题
    )
  
  # 保存为 PDF，文件名包含配色方案的名称
  output_file <- file.path(output_dir, paste0("DimPlot_RColorBrewer_", palette_name, ".pdf"))
  ggsave(
    filename = output_file,
    plot = dimplot,
    width = 8,  # 设置 PDF 宽度
    height = 6  # 设置 PDF 高度
  )
  
  message("已保存: ", output_file)
}

head(seu)

# manual colors ------------------------------------------------------------
names(ggsci_palettes)
palette_function <- ggsci_palettes[["ggsci_npg"]]
# 获取当前配色方案颜色，补足颜色数量
palette_colors <- get_filled_palette(palette_function(n_celltypes), n_celltypes)
palette_colors

max_colors <- brewer.pal.info["Paired", "maxcolors"]
palette_colors <- brewer.pal(n_celltypes, palette_name)
cols <-  c("#E64B35FF", "lightblue", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462" ,"#B3DE69", "#8DD3C7", "#D9D9D9", "#BC80BD", "#CCEBC5")

dimplot <- DimPlot(
  seu,cols= cols, raster = T,pt.size = 5,raster.dpi = c(2048,2048),
  group.by = "celltype",  #  split.by = "geo_accession",
  label = F,      # no显示标签
  label.size = 4      # 标签字体大小
) +
  # scale_color_manual(values = palette_colors) + # 应用配色方案
  theme_minimal() +  # 使用简单的背景主题
  theme(
    legend.position = "right",         # 图例显示在右侧
    panel.grid = element_blank(),      # 去掉网格线
    axis.text = element_blank(),       # 去掉坐标轴文字
    axis.ticks = element_blank(),      # 去掉坐标轴刻度
    axis.title = element_blank(),      # 去掉坐标轴标题
    plot.title = element_blank()       # 去掉 UMAP 图标题
  ); dimplot

ggsave("visual/umap-all.pdf", width = 9, height = 8)

# seu$geo_accession

# splitby=group -----------------------------------------------------------
head(seu)
dimplot <- DimPlot(
  seu,cols= cols, raster = T,
  group.by = "celltype",  split.by = "group",
  label = F     # no显示标签
  # label.size = 4      # 标签字体大小
) +
  # scale_color_manual(values = palette_colors) + # 应用配色方案
  theme_minimal() +  # 使用简单的背景主题
  theme(
    legend.position = "right",         # 图例显示在右侧
    panel.grid = element_blank(),      # 去掉网格线
    axis.text = element_blank(),       # 去掉坐标轴文字
    axis.ticks = element_blank(),      # 去掉坐标轴刻度
    axis.title = element_blank(),      # 去掉坐标轴标题
    plot.title = element_blank()       # 去掉 UMAP 图标题
  )
dimplot
# DimPlot(seu,split.by = "geo_accession", group.by = "celltype")
ggsave("visual/umap-split.pdf", width = 16, height = 8)


# Cell ratio --------------------------------------------------------------
library(ggalluvial)
# 设置分组身份为 "celltype"
Idents(seu) <- "celltype"

# 按照分组列 "group" 计算细胞数和比例
table(Idents(seu), seu$group) # 查看不同细胞类型在每个组中的细胞数
Cellratio <- prop.table(table(Idents(seu), seu$group), margin = 2) # 计算各组不同细胞类型的比例
Cellratio <- as.data.frame(Cellratio)

# 两列变六列 -------------------------------------------------------------------

# 按照分组列 "group" 计算细胞数和比例
table(Idents(seu), seu$orig.ident) # 查看不同细胞类型在每个组中的细胞数
Cellratio <- prop.table(table(Idents(seu), seu$orig.ident), margin = 2) # 计算各组不同细胞类型的比例
Cellratio <- as.data.frame(Cellratio)

# 重命名列名
colnames(Cellratio) <- c("celltype", "group", "ratio")
colnames(Cellratio) <- c("celltype", "sample", "ratio")

# 调整 x 轴分组顺序（可根据实际顺序调整）  已有顺序
# Cellratio$group <- factor(Cellratio$group, levels = c("Group1", "Group2", "Group3")) # 修改成你的组名顺序

# 自定义细胞类型的颜色
# color_cluster <- c("red", "blue", "green", "yellow") # 修改成你需要的颜色

# 绘图
p <- ggplot(Cellratio, aes(x = group, y = ratio, fill = celltype, stratum = celltype, alluvium = celltype)) +
  geom_col(width = 0.4, color = NA) +
  geom_flow(width = 0.4, alpha = 0.2, knot.pos = 0) + # 使用 knot.pos 调整流线的弯曲程度
  scale_fill_manual(values = cols) +
  theme_classic() +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid")) +
  theme(legend.position = 'none') # 取消图例展示
p
# 保存图形
ggsave(plot = p, filename = "visual/celltype_ratio.pdf", width = 4, height = 8)


# 两列变六列 -------------------------------------------------------------------

# 绘图
p <- ggplot(Cellratio, aes(x = sample, y = ratio, fill = celltype, stratum = celltype, alluvium = celltype)) +
  geom_col(width = 0.4, color = NA) +
  geom_flow(width = 0.4, alpha = 0.2, knot.pos = 0) + # 使用 knot.pos 调整流线的弯曲程度
  scale_fill_manual(values = cols) +
  theme_classic() +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid")) +
  theme(legend.position = 'none') # 取消图例展示
p
# 保存图形
ggsave(plot = p, filename = "visual/celltype_ratio.pdf", width = 4, height = 8)

