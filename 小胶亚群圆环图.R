# Load necessary libraries
library(Seurat)
library(dplyr)
library(ggplot2)

# Calculate the proportion of each cell type in each group
proportion_data <- microglia_cells@meta.data %>%
  group_by(seurat_clusters, group) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  group_by(group) %>%
  mutate(proportion = count / sum(count))

# 加载必要的库
library(ggplot2)
library(dplyr)
library(ggsci)  # 使用 ggsci 包中的 npg 配色

# 过滤数据，仅保留 Control 和 AD + WT 两组
filtered_data <- proportion_data
# %>%
#   filter(group %in% c("Control", "AD + WT"))

# 绘制圆环图，使用 ggsci 的 npg 配色
micro_ratio <- ggplot(filtered_data, aes(x = 2, y = proportion, fill = seurat_clusters)) +
  geom_bar(stat = "identity", width = 1, color = "white") +  # 绘制条形图
  coord_polar(theta = "y") +  # 转换为极坐标
  facet_wrap(~ group) +  # 为每个 group 单独生成一个圆环图
  xlim(0.5, 2.5) +  # 设置图形范围，使其成为圆环
  theme_void() +  # 移除背景和坐标轴
  theme(
    legend.position = "right",  # 图例位置
    plot.title = element_text(hjust = 0.5),
    strip.text = element_text(size = 14)  # 调整分面标签的字体大小
  ) +
  # scale_fill_npg() +  # 使用 npg 配色（对比鲜明）
  scale_fill_manual(values = micro_umap_pal) +  # 使用自定义配色方案
  geom_text(aes(label = scales::percent(proportion, accuracy = 0.1)), 
            position = position_stack(vjust = 0.5), size = 4, color = "black") +  # 在环形内显示百分比
  labs(
    fill = "Microglia Subtype",  # 图例标题
    title = "Proportion of Microglia Subtypes by Group"
  )
micro_ratio
ggsave("visual/micro_ratio.pdf",micro_ratio, width = 8, height = 6)

# 批量配色 --------------------------------------------------------------------

# 加载必要的库
library(ggplot2)
library(dplyr)
library(ggsci)

# 确保目标文件夹存在
dir.create("visual/micro_sub", recursive = TRUE, showWarnings = FALSE)

# 过滤数据，仅保留 Control 和 AD + WT 两组
filtered_data <- proportion_data %>%
  filter(group %in% c("Control", "AD + WT"))

# 定义 ggsci 中的配色方案列表
color_schemes <- c("npg", "nejm", "lancet", "jama", "jco", "ucscgb", "d3", "locuszoom", "igv", "uchicago", "aaas")

# 循环生成和保存图表
for (scheme in color_schemes) {
  # 动态选择配色
  plot <- ggplot(filtered_data, aes(x = 2, y = proportion, fill = RNA_snn_res.0.2)) +
    geom_bar(stat = "identity", width = 1, color = "white") +
    coord_polar(theta = "y") +
    facet_wrap(~ group) +
    xlim(0.5, 2.5) +
    theme_void() +
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5),
      strip.text = element_text(size = 14)
    ) +
    geom_text(aes(label = scales::percent(proportion, accuracy = 0.1)), 
              position = position_stack(vjust = 0.5), size = 4, color = "black") +
    labs(
      fill = "Microglia Subtype",
      title = paste("Proportion of Microglia Subtypes by Group -", scheme)
    ) +
    do.call(paste0("scale_fill_", scheme), list())  # 应用当前配色方案
  
  # 保存为 PDF 文件
  ggsave(filename = paste0("visual/micro_sub/", scheme, "_microglia_subtypes.pdf"), 
         plot = plot, width = 8, height = 4)
}



# d3方案 --------------------------------------------------------------------

# 绘制圆环图，使用 d3 配色方案，并仅显示 cluster3 的百分比
plot <- ggplot(filtered_data, aes(x = 2, y = proportion, fill = RNA_snn_res.0.2)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar(theta = "y") +
  facet_wrap(~ group) +
  xlim(0.5, 2.5) +
  theme_void() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5),
    strip.text = element_text(size = 14)
  ) +
  scale_fill_d3() +  # 使用 d3 配色方案
  geom_text(data = subset(filtered_data, RNA_snn_res.0.2 == "3"),  # 仅显示 cluster3 的比例
            aes(label = scales::percent(proportion, accuracy = 0.1)),
            position = position_stack(vjust = 0.5), size = 4, color = "black") +
  labs(
    fill = "Microglia Subtype",
    title = "Proportion of Microglia Subtypes by Group - d3"
  )

# 保存图表为 PDF 文件
ggsave(filename = "visual/micro_sub/d3_microglia_subtypes.pdf", 
       plot = plot, width = 8, height = 4)


# 亚群umap ------------------------------------------------------------------

# 加载必要的库
library(microglia_cellsrat)
library(ggsci)

# 使用 DimPlot 绘制 UMAP，并应用 d3 配色
umap_plot <- DimPlot(microglia_cells, pt.size = 1,
                     reduction = "umap", group.by = "RNA_snn_res.0.2") +
  scale_color_d3() +  # 使用 d3 配色
  labs(title = "UMAP of Microglia Subtypes") +
  theme_minimal(base_line_size = 0)

# 显示 UMAP 图
print(umap_plot)

# 保存 UMAP 图为 PDF
ggsave(filename = "visual/micro_sub/umap_microglia_subtypes_d3.pdf", 
       plot = umap_plot, width = 8, height = 6)

