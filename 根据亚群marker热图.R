# 加载必要的库
library(Seurat)
library(dplyr)
library(pheatmap)
library(viridis)
library(presto)
# 确定每个亚群的top marker基因
marker_list <- list()
for (cluster_id in unique(seu$RNA_snn_res.0.2)) {
  markers <- FindMarkers(seu, ident.1 = cluster_id, group.by = "RNA_snn_res.0.2", only.pos = TRUE)
  top_markers <- markers %>%
    arrange(desc(avg_log2FC)) %>%
    head(3) %>%
    rownames()
  marker_list[[as.character(cluster_id)]] <- top_markers
}

# 提取所有亚群的top markers
all_top_markers <- unique(unlist(marker_list))

# 计算各亚群的top markers的平均表达
avg_exp <- AverageExpression(seu, features = all_top_markers, group.by = "RNA_snn_res.0.2")
avg_exp_mat <- as.matrix(avg_exp$RNA)

# 绘制无聚类的热图
pheatmap(avg_exp_mat,
         color = viridis(100),  # 使用viridis配色
         cluster_rows = FALSE,  # 不进行行聚类
         cluster_cols = FALSE,  # 不进行列聚类
         display_numbers = FALSE,  # 不显示数字
         border_color = "white",  # 增加白色边框以分隔块
         cellwidth = 15, cellheight = 15  # 设置正方形块的尺寸
)



# 增加每个簇的 top markers 数量为 5
marker_list <- list()
for (cluster_id in unique(seu$RNA_snn_res.0.2)) {
  markers <- FindMarkers(seu, ident.1 = cluster_id, group.by = "RNA_snn_res.0.2", only.pos = TRUE)
  top_markers <- markers %>%
    arrange(desc(avg_log2FC)) %>%
    head(5) %>%
    rownames()
  marker_list[[as.character(cluster_id)]] <- top_markers
}

# 提取所有亚群的 top markers
all_top_markers <- unique(unlist(marker_list))

# 计算各亚群的平均表达
avg_exp <- AggregateExpression(seu, features = all_top_markers, group.by = "RNA_snn_res.0.2")
avg_exp_mat <- as.matrix(avg_exp$RNA)

# 绘制无聚类的热图
pheatmap(avg_exp_mat,
         color = viridis(100),  # 使用viridis配色
         cluster_rows = FALSE,  # 不进行行聚类
         cluster_cols = FALSE,  # 不进行列聚类
         display_numbers = FALSE,  # 不显示数字
         border_color = "white",  # 增加白色边框以分隔块
         cellwidth = 15, cellheight = 15  # 设置正方形块的尺寸
)



# 加载必要的库
library(Seurat)

# 提取所有亚群的 top markers
all_top_markers <- unique(unlist(marker_list))

# 使用 DotPlot 显示各亚群的 top marker 基因表达情况
dot_plot <- DotPlot(seu, features = all_top_markers, group.by = "RNA_snn_res.0.2") +
  scale_color_viridis_c() +  # 使用 viridis 配色
  theme_minimal() +
  labs(title = "Expression of Top Markers in Microglia Subtypes")

# 显示 DotPlot
print(dot_plot)



# 加载必要的库
library(Seurat)
library(ggplot2)

# 选择几个代表性的 top marker 基因
top_genes_to_plot <- head(all_top_markers, 9)  # 选择前 9 个基因

# 使用 VlnPlot 显示基因的表达情况
vln_plot <- VlnPlot(seu, features = top_genes_to_plot, group.by = "RNA_snn_res.0.2", 
                    pt.size = 0, combine = TRUE) +
  scale_fill_viridis_d() +  # 使用 viridis 配色
  theme_minimal() +
  labs(title = "Expression of Selected Top Markers in Microglia Subtypes")

# 显示 VlnPlot
print(vln_plot)



# 山鸡图 ---------------------------------------------------------------------

# 创建用于绘图的数据框，包含细胞的分群信息和选择的 top marker 表达值
plot_data <- FetchData(seu, vars = c("RNA_snn_res.0.2", all_top_markers))

# 将数据转换为长格式以适应 ridge plot 的需求
plot_data_long <- plot_data %>%
  pivot_longer(cols = all_top_markers, names_to = "Gene", values_to = "Expression")

# 绘制 Ridge Plot
ridge_plot <- ggplot(plot_data_long, aes(x = Expression, y = RNA_snn_res.0.2, fill = RNA_snn_res.0.2)) +
  geom_density_ridges(scale = 2, rel_min_height = 0.01, color = "white", size = 0.3) +
  scale_fill_viridis_d(option = "D") +  # 使用对比鲜明的 viridis 配色
  facet_wrap(~ Gene, scales = "free_x", ncol = 3) +  # 每个基因独立一列
  theme_minimal() +
  labs(title = "Expression of Top Markers in Microglia Subtypes",
       x = "Expression Level",
       y = "Microglia Subtypes") +
  theme(
    legend.position = "none",  # 移除图例
    strip.text = element_text(size = 10)  # 调整分面标签大小
  )

# 显示 Ridge Plot
print(ridge_plot)



# 重新筛选 --------------------------------------------------------------------

# 加载必要的库
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggridges)

# 重新筛选每个亚群的 top marker 基因，增加表达量的限制条件
marker_list <- list()
for (cluster_id in unique(seu$RNA_snn_res.0.2)) {
  markers <- FindMarkers(seu, ident.1 = cluster_id, group.by = "RNA_snn_res.0.2", only.pos = TRUE)
  top_markers <- markers %>%
    filter(avg_log2FC > 0.5 & pct.1 > 0.2) %>%  # 增加筛选条件
    arrange(desc(avg_log2FC)) %>%
    head(3) %>%
    rownames()
  marker_list[[as.character(cluster_id)]] <- top_markers
}

# 将所有 top markers 合并为一个列表
all_top_markers <- unique(unlist(marker_list))

# 提取每个细胞的分群信息和选择的 top marker 表达值
plot_data <- FetchData(seu, vars = c("RNA_snn_res.0.2", all_top_markers))

# 将数据转换为长格式
plot_data_long <- plot_data %>%
  pivot_longer(cols = all_top_markers, names_to = "Gene", values_to = "Expression")

# 绘制 Ridge Plot
ridge_plot <- ggplot(plot_data_long, aes(x = Expression, y = RNA_snn_res.0.2, fill = RNA_snn_res.0.2)) +
  geom_density_ridges(scale = 2, rel_min_height = 0.01, color = "white", size = 0.3) +
  scale_fill_viridis_d(option = "D") +  # 使用对比鲜明的 viridis 配色
  facet_wrap(~ Gene, scales = "free_x", ncol = 3) +  # 每个基因独立一列
  theme_minimal() +
  labs(title = "Expression of Top Markers in Microglia Subtypes",
       x = "Expression Level",
       y = "Microglia Subtypes") +
  theme(
    legend.position = "none",  # 移除图例
    strip.text = element_text(size = 10)  # 调整分面标签大小
  )

# 显示 Ridge Plot
print(ridge_plot)

