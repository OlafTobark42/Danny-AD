# 加载必要的库
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggridges)
table(seu$celltype)
# 设置分组标识符为 RNA_snn_res.0.2 列
Idents(seu) <- "RNA_snn_res.0.2"

# 使用 FindAllMarkers 查找所有亚群的标记基因
all_markers <- FindAllMarkers(seu, only.pos = TRUE)

# 为每个亚群选择 top 3 个标记基因
top_markers <- all_markers %>%
  group_by(cluster) %>%
  filter(avg_log2FC > 0.5 & pct.1 > 0.2) %>%  # 增加筛选条件
  top_n(3, avg_log2FC) %>%
  ungroup() %>%
  pull(gene) %>%
  unique()

# 加载必要的库
library(Seurat)
library(ggplot2)

# 绘制18个top marker基因的Feature Plot
feature_plot <- FeaturePlot(seu, features = top_markers, cols = c("lightgrey", "blue"), ncol = 3)

# 显示Feature Plot
print(feature_plot)
ggsave("visual/sub_top_marker.pdf", feature_plot,
       width = 12, height = 20)



# 提取所有选择的 top marker 的表达数据
plot_data <- FetchData(seu, vars = c("RNA_snn_res.0.2", top_markers))

# 将数据转换为长格式以适应 ridge plot 的需求
plot_data_long <- plot_data %>%
  pivot_longer(cols = top_markers, names_to = "Gene", values_to = "Expression")

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


# heatmap -----------------------------------------------------------------
# 加载必要的库
library(Seurat)
library(dplyr)
library(pheatmap)
library(Matrix)

# 设置 Idents 为 RNA_snn_res.0.2
Idents(seu) <- "RNA_snn_res.0.2"

# 使用 AverageExpression 提取 top_markers 基因的平均表达
avg_exp <- AverageExpression(seu, features = top_markers, group.by = "RNA_snn_res.0.2")
avg_exp_mat <- as.matrix(avg_exp$RNA)

# 对表达矩阵进行 log1p 变换，以减少极端高值的影响
avg_exp_mat <- log1p(avg_exp_mat)

# 确定每个基因的高表达亚群（找到表达最高的列）
max_cluster <- apply(avg_exp_mat, 1, function(x) colnames(avg_exp_mat)[which.max(x)])

# 将基因和高表达亚群信息组合成数据框
gene_order_df <- data.frame(Gene = rownames(avg_exp_mat), MaxCluster = max_cluster, stringsAsFactors = FALSE)

# 按高表达亚群分组，并在每个亚群内按表达量降序排序
gene_order_df <- gene_order_df %>%
  group_by(MaxCluster) %>%
  arrange(-apply(avg_exp_mat, 1, max)) %>%
  ungroup()

# 根据新的基因顺序重新排列表达矩阵
avg_exp_mat <- avg_exp_mat[gene_order_df$Gene, ]

# 对每个基因（行）进行 Z-score 标准化
avg_exp_mat <- t(scale(t(avg_exp_mat)))

# 绘制排序后的热图
pheatmap(avg_exp_mat,
         color = viridis::viridis(100),
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         display_numbers = FALSE,
         border_color = "white",
         cellwidth = 15, cellheight = 15,
         main = "Ordered Expression of Top Markers in Microglia Subtypes",
         angle_col = 45)


# 指定基因 --------------------------------------------------------------------

# 加载必要的库
library(Seurat)
library(pheatmap)
library(Matrix)

# 设置 Idents 为 RNA_snn_res.0.2
Idents(seu) <- "RNA_snn_res.0.2"

# 指定基因集
selected_genes <- c("P2RY12", "P2RY13", "CX3CR1", "TMEM119", "CD9", "ITGAX", "CLEC7A", "CD63", 
                    "SPP1", "LPL", "TREM2", "APOE", "NAMPT", "ACSL1", "DPYD", "CD163", 
                    "HOMEO", "DAM", "LDAM")
selected_genes <- stringr::str_to_title(selected_genes)

# 将 selected_genes 设为因子，并确保顺序保持一致
selected_genes <- factor(selected_genes, levels = selected_genes)

# 使用 AverageExpression 提取基因的平均表达
avg_exp <- AverageExpression(seu, features = selected_genes, group.by = "RNA_snn_res.0.2")
avg_exp_mat <- as.matrix(avg_exp$RNA)

# 对表达矩阵进行 log1p 变换
avg_exp_mat <- log1p(avg_exp_mat)

# 确保 avg_exp_mat 的行顺序与 selected_genes 一致
avg_exp_mat <- avg_exp_mat[as.character(selected_genes), ]

# 绘制热图，使用 viridis 配色以突出特异性高表达
# 绘制热图，显式关闭行列聚类
pheatmap(avg_exp_mat,
         color = viridis::viridis(100),
         cluster_rows = FALSE,  # 关闭行聚类
         cluster_cols = FALSE,  # 关闭列聚类
         display_numbers = FALSE,
         border_color = "white",
         cellwidth = 15, cellheight = 15,
         main = "Expression of Selected Genes in Microglia Subtypes",
         angle_col = 45)


# 过滤gene ------------------------------------------------------------------
# 加载必要的库
library(Seurat)
library(pheatmap)
library(Matrix)
library(stringr)

# 指定基因集并将基因名称转换为首字母大写格式
selected_genes <- c("P2RY12", "P2RY13", "CX3CR1", "TMEM119", "CD9", "ITGAX", "CLEC7A", "CD63", 
                    "SPP1", "LPL", "TREM2", "APOE", "NAMPT", "ACSL1", "DPYD", "CD163", 
                    "HOMEO", "DAM", "LDAM")
selected_genes <- str_to_title(selected_genes)

# 确保 selected_genes 是因子，且顺序保持一致
selected_genes <- factor(selected_genes, levels = selected_genes)

# 检查基因是否在 Seurat 对象中存在
existing_genes <- selected_genes[selected_genes %in% rownames(seu@assays$RNA$data)]
existing_genes
# 计算选择基因在各亚群的平均表达
Idents(seu) <- "RNA_snn_res.0.2"
avg_exp <- AverageExpression(seu, features = existing_genes, group.by = "RNA_snn_res.0.2")
avg_exp_mat <- as.matrix(avg_exp$RNA)

# 对表达矩阵进行 log1p 变换，以减少极端高值的影响
avg_exp_mat <- log1p(avg_exp_mat)

# 确保 avg_exp_mat 的行顺序与 selected_genes 一致
avg_exp_mat <- avg_exp_mat[as.character(existing_genes), ]

# 对每个基因（行）进行 Z-score 标准化
avg_exp_mat <- t(scale(t(avg_exp_mat)))

# 绘制横向热图，调整字体大小和单元格宽高
ht <- pheatmap(avg_exp_mat,
               color = viridis::viridis(100),
               cluster_rows = FALSE,  # 关闭行聚类
               cluster_cols = FALSE,  # 关闭列聚类
               display_numbers = FALSE,
               border_color = "white",
               cellwidth = 15, cellheight = 15,  # 控制单元格尺寸
               fontsize_row = 8,  # 控制基因名字体大小
               fontsize_col = 10, # 控制亚群名字体大小
               main = "Expression of Selected Genes in Microglia Subtypes",
               angle_col = 45)  # 倾斜列名
print(ht)
ggsave("visual/sub_heat_fromNature_shu.pdf",ht,
       width = 12, height = 8)




 # 另一种基因集 ------------------------------------------------------------------

# 加载必要的库
library(Seurat)
library(pheatmap)
library(stringr)

# 指定新的基因集并将基因名称转换为首字母大写格式
selected_genes <- c("Cacnb2", "Tmem119", "Zfhx3", "P2ry12", "Maf", "Rtp4", "Stat1", "Trim30a", 
                    "Ifitm3", "Ccl12", "Cxcl10", "Ifit3", "Isg15", "Apoe", "Cacna1a", "Cst7", 
                    "Etl4", "Ccl4", "Ccl3", "Myo1e", "Lpl", "Cd14", "Ctnna3", "Lrrtm3", "Ank", 
                    "Mamdc2", "Wfdc17", "Ftl1", "Fth1", "Rps24", "Spp1", "Gpnmb", "H2-Aa", 
                    "Cd74", "H2-Ab1", "H2-Eb1", "H2-Q7")
selected_genes <- str_to_title(selected_genes)

# 仅提取在 Seurat 对象中存在的基因
existing_genes <- selected_genes[selected_genes %in% rownames(seu@assays$RNA$data)]

# 使用 FetchData 获取表达数据
Idents(seu) <- "RNA_snn_res.0.2"
avg_exp <- AverageExpression(seu, features = existing_genes, group.by = "RNA_snn_res.0.2")
avg_exp_mat <- as.matrix(avg_exp$RNA)

# 对表达矩阵进行 log1p 变换，以减少极端高值的影响
avg_exp_mat <- log1p(avg_exp_mat)

# 确保 avg_exp_mat 的行顺序与 selected_genes 一致
avg_exp_mat <- avg_exp_mat[as.character(existing_genes), ]

# 对每个基因（行）进行 Z-score 标准化
avg_exp_mat <- t(scale(t(avg_exp_mat)))

# 转置矩阵以便基因名显示在 x 轴
avg_exp_mat <- t(avg_exp_mat)

# 绘制横向热图，调整字体大小和单元格宽高
ht <- pheatmap(avg_exp_mat,
         color = viridis::viridis(100),
         cluster_rows = FALSE,  # 关闭行聚类
         cluster_cols = FALSE,  # 关闭列聚类
         display_numbers = FALSE,
         border_color = "white",
         cellwidth = 15, cellheight = 15,  # 控制单元格尺寸
         fontsize_row = 8,  # 控制基因名字体大小
         fontsize_col = 10, # 控制亚群名字体大小
         main = "Expression of Selected Genes in Microglia Subtypes",
         angle_col = 45)  # 倾斜列名
print(ht)
ggsave("visual/sub_heat.pdf",ht,
       width = 12, height = 8)


