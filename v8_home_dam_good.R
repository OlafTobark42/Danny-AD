rm(list = ls())
library(ggplot2)
library(viridis)
# Load necessary libraries
library(Seurat)
library(pheatmap)
library(stringr)

microglia_cells <- readRDS("rds/micro_sub.rds")
# UMAP图展示小胶质细胞亚群
DimPlot(microglia_cells, reduction = "umap", label = TRUE) + ggtitle("Microglia Subclusters")

library(ggrastr); library(Nebulosa)
ggrastr::rasterize(Nebulosa::plot_density(microglia_cells,
                                          c("Mertk","Pparg"),
                                          size = 0.2), dpi = 300)
FeaturePlot(microglia_cells, c("Mertk","Pparg"))

# 按样本 `orig.ident` 分组的 UMAP 图
DimPlot(microglia_cells, reduction = "umap", ncol = 3,
        split.by = "group", label = T) + ggtitle("UMAP of Microglia Subclusters by Group")


# 找出每个小胶质细胞亚群的特异性marker
# microglia_markers <- FindAllMarkers(microglia_cells, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)


# List of marker genes (example: from Homeostatic, DAM, LDAM categories)
marker_genes <- str_to_title(c("P2RY12", "P2RY13", "CX3CR1", "TMEM119", 
                               "CD9", "ITGAX", "CLEC7A", "CD63", "SPP1", "LPL", "TREM2", "APOE",
                               "NAMPT", "ACSL1", "DPYD", "CD163"))

home_marker <- c("P2ry12", "P2ry13", "Cx3cr1", "Tmem119")
DAM_marker <- c("Cd9", "Itgax", "Cd63", "Spp1", "Lpl", "Trem2", "Apoe")
LDAM_marker <- c("Nampt", "Acsl1", "Dpyd", "Cd163")
LDAM_marker %in% rownames(microglia_cells)

microglia_cells <- AddModuleScore(microglia_cells, 
                                  features = list(DAM_marker), name = "DAM")
microglia_cells <- AddModuleScore(microglia_cells, 
                                  features = list(LDAM_marker), name = "LDAM")
microglia_cells <- AddModuleScore(microglia_cells, 
                                  features = list(home_marker), name = "Home")
head(microglia_cells)

# Home1
mydata <- FetchData(microglia_cells, vars = c("umap_1", "umap_2", "Home1"))

# 计算颜色梯度的上下限
Home1_min <- quantile(mydata$Home1, 0.01)  # 设置5%作为最低点
Home1_max <- quantile(mydata$Home1, 0.95)  # 设置95%作为最高点

# 调整颜色梯度的上下限并添加梯度线阈值
ggplot(mydata, aes(x = umap_1, y = umap_2, colour = Home1)) +
  geom_point(size = 1) +  # 调整点的大小
  scale_color_viridis(option = "viridis") +  # 限制颜色梯度会导致灰色点太多 , limits = c(Home1_min, Home1_max)
  theme_minimal() +  # 使用简约主题
  theme(
    legend.position = "none",  # 完全移除 legend
    panel.background = element_rect(fill = "black"),  # 设置黑色背景
    panel.grid.major = element_blank(),  # 移除网格线
    panel.grid.minor = element_blank(),  # 移除小网格线
    plot.background = element_rect(fill = "black"),  # 设置黑色绘图区背景
    axis.title = element_text(color = "white"),  # 坐标轴标题颜色
    axis.text = element_blank()  # 移除坐标轴刻度
  ) +
  geom_density_2d(data = subset(mydata, Home1 > Home1_max * 0.6), color = "white")   # 添加白色梯度线（仅在高亮区域）



# DAM
mydata <- FetchData(microglia_cells, vars = c("umap_1", "umap_2", "DAM1"))

# 计算颜色梯度的上下限
Home1_min <- quantile(mydata$DAM1, 0.01)  # 设置5%作为最低点
Home1_max <- quantile(mydata$DAM1, 0.95)  # 设置95%作为最高点

# 调整颜色梯度的上下限并添加梯度线阈值
ggplot(mydata, aes(x = umap_1, y = umap_2, colour = DAM1)) +
  geom_point(size = 1) +  # 调整点的大小
  scale_color_viridis(option = "viridis") +  # 限制颜色梯度会导致灰色点太多 , limits = c(Home1_min, Home1_max)
  theme_minimal() +  # 使用简约主题
  theme(
    legend.position = "none",  # 完全移除 legend
    panel.background = element_rect(fill = "black"),  # 设置黑色背景
    panel.grid.major = element_blank(),  # 移除网格线
    panel.grid.minor = element_blank(),  # 移除小网格线
    plot.background = element_rect(fill = "black"),  # 设置黑色绘图区背景
    axis.title = element_text(color = "white"),  # 坐标轴标题颜色
    axis.text = element_blank()  # 移除坐标轴刻度
  ) +
  geom_density_2d(data = subset(mydata, DAM1 > Home1_max * 0.85), color = "white")   # 添加白色梯度线（仅在高亮区域）
# 0.85 better than 0.7 and 0.9


# LDAM
mydata <- FetchData(microglia_cells, vars = c("umap_1", "umap_2", "LDAM1"))

# 计算颜色梯度的上下限
Home1_min <- quantile(mydata$LDAM1, 0.01)  # 设置5%作为最低点
Home1_max <- quantile(mydata$LDAM1, 0.95)  # 设置95%作为最高点

# 调整颜色梯度的上下限并添加梯度线阈值
ggplot(mydata, aes(x = umap_1, y = umap_2, colour = LDAM1)) +
  geom_point(size = 1) +  # 调整点的大小
  scale_color_viridis(option = "viridis") +  # 限制颜色梯度会导致灰色点太多 , limits = c(Home1_min, Home1_max)
  theme_minimal() +  # 使用简约主题
  theme(
    legend.position = "none",  # 完全移除 legend
    panel.background = element_rect(fill = "black"),  # 设置黑色背景
    panel.grid.major = element_blank(),  # 移除网格线
    panel.grid.minor = element_blank(),  # 移除小网格线
    plot.background = element_rect(fill = "black"),  # 设置黑色绘图区背景
    axis.title = element_text(color = "white"),  # 坐标轴标题颜色
    axis.text = element_blank()  # 移除坐标轴刻度
  ) +
  geom_density_2d(data = subset(mydata, LDAM1 > Home1_max * 0.9), color = "white")



# manual add subtype
microglia_cells$subtype <- NA
microglia_cells$subtype <- ifelse(microglia_cells$seurat_clusters %in% 3, "DAM", "HOMEOSTATIC")

table(microglia_cells$subtype, microglia_cells$group)

# 加载必要的库
library(ggplot2)
library(dplyr)

# 构建数据框并设置分组顺序
data <- data.frame(
  subtype = rep(c("DAM", "HOMEOSTATIC"), each = 3),
  group = rep(c("Control", "AD + WT", "AD + Mertk-/-"), times = 2),
  count = c(64, 392, 808, 7582, 5035, 6698)
)

data$group <- factor(data$group, levels = c("Control", "AD + WT", "AD + Mertk-/-"))

# 计算每个 group 的总细胞数
data <- data %>%
  group_by(group) %>%
  mutate(total = sum(count)) %>%
  ungroup() %>%
  mutate(proportion = count / total)

# 自定义颜色方案
custom_colors <- c("DAM" = "#66c2a5", "HOMEOSTATIC" = "#fc8d62")

# 绘制圆环图
ggplot(data, aes(x = 2, y = proportion, fill = subtype)) +
  geom_bar(stat = "identity", width = 1, color = "white") +  # 绘制条形图
  coord_polar(theta = "y") +  # 将条形图转换为极坐标
  facet_wrap(~ group) +  # 每个 group 单独一张图
  xlim(0.5, 2.5) +  # 调整图形的半径范围，使其看起来像环形
  theme_void() +  # 移除背景和坐标轴
  theme(
    legend.position = "right",  # 图例位置
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(
    fill = "Subtype",  # 图例标题
    title = "Proportion of Microglia Subtypes by Group"
  ) +
  scale_fill_manual(values = custom_colors)  # 应用自定义颜色

# 修改配色
# 构建数据框并设置分组顺序
data <- data.frame(
  subtype = rep(c("DAM", "HOMEOSTATIC"), each = 3),
  group = rep(c("Control", "AD + WT", "AD + Mertk-/-"), times = 2),
  count = c(64, 392, 808, 7582, 5035, 6698)
)

# 将 subtype 设为 factor，并确保 HOMEOSTATIC 在前
data$subtype <- factor(data$subtype, levels = c("HOMEOSTATIC", "DAM"))
data$group <- factor(data$group, levels = c("Control", "AD + WT", "AD + Mertk-/-"))

# 计算每个 group 的总细胞数
data <- data %>%
  group_by(group) %>%
  mutate(total = sum(count)) %>%
  ungroup() %>%
  mutate(proportion = count / total)

# 自定义颜色方案（灰色和黄色）
custom_colors <- c("HOMEOSTATIC" = "#bdbdbd", "DAM" = "#fdae61")

# 绘制圆环图
ggplot(data, aes(x = 2, y = proportion, fill = subtype)) +
  geom_bar(stat = "identity", width = 1, color = "white") +  # 绘制条形图
  coord_polar(theta = "y") +  # 将条形图转换为极坐标
  facet_wrap(~ group) +  # 每个 group 单独一张图
  xlim(0.5, 2.5) +  # 调整图形的半径范围，使其看起来像环形
  theme_void() +  # 移除背景和坐标轴
  theme(
    legend.position = "right",  # 图例位置
    plot.title = element_text(hjust = 0.5),
    strip.text = element_text(size = 14)  # 调整分面标签的字体大小
  ) +
  scale_fill_manual(values = custom_colors) +  # 应用自定义颜色
  geom_text(data = subset(data, subtype == "DAM"), 
            aes(label = scales::percent(proportion, accuracy = 0.1)),
            position = position_stack(vjust = 0.5), size = 4, color = "black") +  # 在圆环中心添加比例
  labs(
    fill = "Subtype",  # 图例标题
    title = "Proportion of Microglia Subtypes by Group"
  )


