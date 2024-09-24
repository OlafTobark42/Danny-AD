rm(list = ls())
library(Seurat)

# ReadRDS
seurat_object = readRDS("rds/micro.rds")

# 定义基因集

homeostatic_genes <- c("P2ry12", "P2ry13", "Cx3cr1", "Tmem119")
DAM_genes <- c("Cd9", "Itgax", "Clec7a", "Cd63", "Spp1", "Lpl", "Trem2", "Apoe")
LDAM_genes <- c("Nampt", "Acsl1", "Dpyd", "Cd163")

# 将基因集组合成列表
gene_sets <- list(
  Homeostatic = homeostatic_genes,
  DAM = DAM_genes,
  LDAM = LDAM_genes
)

# 假设Seurat对象名为seurat_object
# 使用AddModuleScore打分
seurat_object <- AddModuleScore(
  object = seurat_object,
  features = gene_sets,
  name = c("Homeostatic_Score", "DAM_Score", "LDAM_Score")
)

# 可视化Homeostatic分数
FeaturePlot(seurat_object, features = "Homeostatic_Score1", cols = c("lightgrey", "blue"))

# 可视化DAM分数
FeaturePlot(seurat_object, features = "DAM_Score1", cols = c("lightgrey", "red"))

# 可视化LDAM分数
FeaturePlot(seurat_object, features = "LDAM_Score1", cols = c("lightgrey", "green"))

# 查看各cluster的分数分布，使用小提琴图
VlnPlot(seurat_object, features = c("Homeostatic_Score1", "DAM_Score1", "LDAM_Score1"), group.by = "seurat_clusters")

# 手动注释cluster为亚群
seurat_object$cell_type <- Idents(seurat_object)
seurat_object$cell_type <- factor(seurat_object$cell_type, levels = c("Cluster1", "Cluster2", "Cluster3", "Cluster4", "Cluster5"))

# 假设你根据分数注释：
seurat_object$cell_type[seurat_object$seurat_clusters == "1"] <- "Homeostatic Microglia"
seurat_object$cell_type[seurat_object$seurat_clusters == "2"] <- "DAM Microglia"
seurat_object$cell_type[seurat_object$seurat_clusters == "3"] <- "LDAM Microglia"
# 为其余的cluster注释

# 查看细胞类型的UMAP分布
DimPlot(seurat_object, group.by = "cell_type", label = TRUE)

# 加载必要的包
library(Seurat)
library(ggplot2)

# 假设 Seurat 对象是 seurat_object
# 设置黑底的主题
umap_black <- DimPlot(seurat_object, reduction = "umap", group.by = "cell_type") +
  theme(panel.background = element_rect(fill = "black", color = "black"),
        plot.background = element_rect(fill = "black", color = "black"),
        axis.line = element_line(color = "white"),
        axis.ticks = element_line(color = "white"),
        axis.text = element_text(color = "white"),
        legend.text = element_text(color = "white"),
        legend.title = element_text(color = "white"),
        panel.grid = element_blank()) +
  scale_color_viridis_d(option = "plasma")  # 选择合适的颜色方案

# 显示黑底 UMAP
print(umap_black)

# 安装和加载pheatmap或使用DoHeatmap
install.packages("pheatmap")
library(pheatmap)

# 提取用于绘制热图的基因表达矩阵
# 假设你已经定义了三个基因集，如之前提到的Homeostatic, DAM, LDAM
gene_list <- c(homeostatic_genes, DAM_genes, LDAM_genes)

# 从Seurat对象中提取表达矩阵
expression_matrix <- GetAssayData(seurat_object, slot = "scale.data")[gene_list, ]

# 绘制热图 (使用virdis颜色)
pheatmap(expression_matrix,
         color = viridis::viridis(100),  # 蓝-绿配色
         cluster_rows = TRUE,  # 根据基因分群
         cluster_cols = TRUE,  # 根据细胞分群
         scale = "row",  # 对行进行z-score归一化
         show_rownames = TRUE,  # 显示基因名字
         show_colnames = FALSE)  # 不显示细胞名字

