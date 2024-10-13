seu <- readRDS("rds/manual_annoted.rds")
head(seu);table(seu$celltype)
DimPlot(seu, group.by = "celltype", label = T)
DimPlot(seu, group.by = "seurat_clusters", label = T)

Idents(seu) <- seu$celltype
deg_20 <- FindMarkers(seu, ident.1 = 20, group.by = "seurat_clusters")
deg_21 <- FindMarkers(seu, ident.1 = 21, group.by = "seurat_clusters")
deg_22 <- FindMarkers(seu, ident.1 = 22, group.by = "seurat_clusters")
write.csv(deg_22,"deg_result/unknwon_22.csv")

seu$celltype <- ifelse(seu$seurat_clusters %in% 25, "Macrophages", seu$celltype)
seu$celltype <- ifelse(seu$seurat_clusters %in% 20, "Immune cells", seu$celltype)
seu$celltype <- ifelse(seu$seurat_clusters %in% 22, "Astrocytes", seu$celltype)
seu$celltype <- ifelse(seu$seurat_clusters %in% 21, "Fibroblasts", seu$celltype)


# Load the necessary libraries
library(Seurat)
library(ggplot2)
library(ggsci)  # for ggsci color palettes

# Reorder the cell types: Microglia first, Neuron last
seu$celltype <- factor(seu$celltype, levels = c("Microglia", "Oligodendrocytes", "Astrocytes", 
                                                "Endothelial Cells", "Ependymal Cells", "VSMCs", 
                                                "Macrophages", "Fibroblasts", "Immune cells", "Neuron"))
saveRDS(seu, "rds/done_annotated.rds")
table(seu$celltype, seu$group)
# 加载必要的库
library(ggplot2)
library(dplyr)

# 构建数据框
data <- data.frame(
  celltype = rep(c("Microglia", "Oligodendrocytes", "Astrocytes", "Endothelial Cells", 
                   "Ependymal Cells", "VSMCs", "Macrophages", "Fibroblasts", 
                   "Immune cells", "Neuron"), each = 3),
  group = rep(c("Control", "AD + WT", "AD + Mertk-/-"), times = 10),
  count = c(7646, 5427, 7506, 8761, 4926, 9319, 1465, 1209, 977, 5650, 3464, 4172,
            726, 250, 300, 338, 252, 259, 808, 567, 727, 140, 85, 93, 196, 75, 74, 776, 885, 803)
)

# 计算每个 group 的总细胞数
data <- data %>%
  group_by(group) %>%
  mutate(total = sum(count)) %>%
  ungroup() %>%
  mutate(proportion = count / total)

# 自定义颜色
custom_colors <- c("Microglia" = "#d73027", "Oligodendrocytes" = "#fc8d59", "Astrocytes" = "#fee090",
                   "Endothelial Cells" = "#91bfdb", "Ependymal Cells" = "#4575b4", "VSMCs" = "#313695",
                   "Macrophages" = "#a6d96a", "Fibroblasts" = "#66bd63", "Immune cells" = "#1a9850", 
                   "Neuron" = "#f46d43")

# 绘制堆叠柱状图
ggplot(data, aes(x = group, y = proportion, fill = celltype)) +
  geom_bar(stat = "identity", color = "white", width = 0.8) +  # 堆叠柱状图
  scale_fill_manual(values = custom_colors) +  # 应用自定义颜色
  theme_minimal() +  # 使用简约主题
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # 旋转X轴标签
    legend.position = "right",  # 图例位置
    legend.title = element_text(size = 12),  # 图例标题大小
    legend.text = element_text(size = 10),  # 图例文字大小
    panel.grid.major = element_blank(),  # 移除主网格线
    panel.grid.minor = element_blank(),  # 移除次网格线
    axis.title.x = element_blank()  # 移除X轴标题
  ) +
  labs(
    y = "Cell Proportion",  # Y轴标签
    fill = "Cell Type"  # 图例标题
  )

# g g s ci nog
# 加载必要的库
library(ggplot2)
library(dplyr)
library(ggsci)  # 用于npg配色

# 构建数据框
data <- data.frame(
  celltype = rep(c("Microglia", "Oligodendrocytes", "Astrocytes", "Endothelial Cells", 
                   "Ependymal Cells", "VSMCs", "Macrophages", "Fibroblasts", 
                   "Immune cells", "Neuron"), each = 3),
  group = rep(c("Control", "AD + WT", "AD + Mertk-/-"), times = 10),
  count = c(7646, 5427, 7506, 8761, 4926, 9319, 1465, 1209, 977, 5650, 3464, 4172,
            726, 250, 300, 338, 252, 259, 808, 567, 727, 140, 85, 93, 196, 75, 74, 776, 885, 803)
)

# 将 celltype 设为 factor，并定义顺序
data$celltype <- factor(data$celltype, levels = c("Microglia", "Oligodendrocytes", "Astrocytes", 
                                                  "Endothelial Cells", "Ependymal Cells", 
                                                  "VSMCs", "Macrophages", "Fibroblasts", 
                                                  "Immune cells", "Neuron"))

data$group <- factor(data$group, levels = c("Control", "AD + WT", "AD + Mertk-/-"))

# 计算每个 group 的总细胞数
data <- data %>%
  group_by(group) %>%
  mutate(total = sum(count)) %>%
  ungroup() %>%
  mutate(proportion = count / total)

# 绘制堆叠柱状图，使用ggsci的npg配色
ggplot(data, aes(x = group, y = proportion, fill = celltype)) +
  geom_bar(stat = "identity", color = "white", width = 0.8) +  # 堆叠柱状图
  scale_fill_npg() +  # 使用npg配色
  theme_minimal() +  # 使用简约主题
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # 旋转X轴标签
    legend.position = "right",  # 图例位置
    legend.title = element_text(size = 12),  # 图例标题大小
    legend.text = element_text(size = 10),  # 图例文字大小
    panel.grid.major = element_blank(),  # 移除主网格线
    panel.grid.minor = element_blank(),  # 移除次网格线
    axis.title.x = element_blank()  # 移除X轴标题
  ) +
  labs(
    y = "Cell Proportion",  # Y轴标签
    fill = "Cell Type"  # 图例标题
  )

