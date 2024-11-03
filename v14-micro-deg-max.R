library(Seurat)
library(dplyr)
library(ggplot2)

library(Seurat)
library(dplyr)
library(ggplot2)

# 获取所有细胞类型
cell_types <- unique(seu$celltype)

# 初始化一个数据框用于存储DEG数量
deg_counts <- data.frame(celltype = character(), DEG_count = integer())
Idents(seu) <- seu$celltype
# 循环计算每个细胞类型的DEG数量
for (cell_type in cell_types) {
  # Subset Seurat object for the specific cell type
  subset_seu <- subset(seu, idents = cell_type)
  Idents(subset_seu) <- subset_seu$group
  # 计算差异表达基因
  deg_result <- FindMarkers(subset_seu, ident.1 = "AD + WT", ident.2 = "Control")
  
  # 筛选出符合条件的DEG (p_val < 0.05 且 avg_log2FC > 1 或 avg_log2FC < -1)
  filtered_deg <- deg_result %>%
    filter(p_val < 0.05 & (avg_log2FC > 1 | avg_log2FC < -1))
  
  # 记录筛选后DEG的数量
  deg_count <- nrow(filtered_deg)
  deg_counts <- rbind(deg_counts, data.frame(celltype = cell_type, DEG_count = deg_count))
}

# 绘制条形图
ggplot(deg_counts, aes(x = reorder(celltype, -DEG_count), y = DEG_count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  labs(title = "Number of Significant DEGs per Cell Type (AD + WT vs Control)", 
       x = "Cell Type", 
       y = "DEG Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



# 多重火山图 -------------------------------------------------------------------

library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(ggsci)
library(ggrepel)

# 获取所有细胞类型
cell_types <- unique(seu$celltype)

# 创建一个空列表来存储每个细胞类型的火山图
volcano_plots <- list()

# 循环遍历每个细胞类型并绘制火山图
for (cell_type in cell_types) {
  # Subset Seurat object for the specific cell type
  subset_seu <- subset(seu, idents = cell_type)
  Idents(subset_seu) <- subset_seu$group
  # 计算差异表达基因
  deg_result <- FindMarkers(subset_seu, ident.1 = "AD + WT", ident.2 = "Control")
  
  # 筛选出符合条件的DEG (adjusted p_val < 0.05)
  deg_result <- deg_result %>%
    filter(p_val_adj < 0.05)
  deg_result$gene <- rownames(deg_result)
  # 创建火山图
  volcano_plot <- ggplot(deg_result, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
    geom_point(aes(color = p_val_adj < 0.01), size = 1.5) +
    scale_color_manual(values = c("black", "red")) +
    geom_text_repel(data = deg_result %>%
                      filter(abs(avg_log2FC) > 1.5 & p_val_adj < 0.01),
                    aes(label = gene), size = 3, max.overlaps = 10) +
    theme_minimal() +
    labs(title = paste("Volcano Plot -", cell_type),
         x = "Log2 Fold Change", 
         y = "-Log10 Adjusted P-Value") +
    theme(legend.position = "none")
  
  # 将火山图添加到列表中
  volcano_plots[[cell_type]] <- volcano_plot
}

# 使用 cowplot 将所有火山图组合到一起
combined_plot <- plot_grid(plotlist = volcano_plots, ncol = 3)

# 显示组合后的火山图
print(combined_plot)

# 修改scRNAvol包函数 -----------------------------------------------------------

library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(purrr)
library(cowplot)

# 假设seu是你的Seurat对象，已经包含了所有细胞类型和样本组别
# 获取所有细胞类型
cell_types <- unique(seu$celltype)

# 初始化一个空的数据框，用于存储所有细胞类型的DEG
combined_deg_result <- data.frame()

# 循环计算每个细胞类型的DEG并合并
for (cell_type in cell_types) {
  # Subset Seurat object for the specific cell type
  subset_seu <- subset(seu, idents = cell_type)
  Idents(subset_seu) <- subset_seu$group
  # 计算差异表达基因
  deg_result <- FindMarkers(subset_seu, ident.1 = "AD + WT", ident.2 = "Control")
  
  # 添加细胞类型信息
  deg_result$cluster <- cell_type
  deg_result$gene <- rownames(deg_result)
  # 筛选出符合条件的DEG (p_val < 0.05 且 avg_log2FC > 1 或 avg_log2FC < -1)
  filtered_deg <- deg_result %>%
    filter(p_val < 0.05 & (avg_log2FC > 1 | avg_log2FC < -1))
  
  # 合并到总数据框中
  combined_deg_result <- rbind(combined_deg_result, filtered_deg)
}

# 调整jjVolcano函数以适应合并后的DEG数据
jjVolcano_combined <- function(
    diffData = combined_deg_result,
    log2FC.cutoff = 0.25,
    pvalue.cutoff = 0.05,
    adjustP.cutoff = 0.01,
    topGeneN = 5,
    col.type = "updown",
    aesCol = c("#0099CC", "#CC3333"),
    legend.position = c(0.7, 0.9),
    base_size = 14,
    tile.col = scales::hue_pal()(length(unique(diffData$cluster))),
    cluster.order = NULL,
    ...) {
  
  # 过滤数据
  diff.marker <- diffData %>%
    filter(abs(avg_log2FC) >= log2FC.cutoff & p_val < pvalue.cutoff)
  
  # 分配类型
  diff.marker <- diff.marker %>%
    mutate(type = ifelse(avg_log2FC >= log2FC.cutoff, "sigUp", "sigDown")) %>%
    mutate(type2 = ifelse(p_val_adj < adjustP.cutoff,
                          paste("adjust Pvalue < ", adjustP.cutoff, sep = ""),
                          paste("adjust Pvalue >= ", adjustP.cutoff, sep = "")
    ))
  
  # 确定cluster的顺序
  if (!is.null(cluster.order)) {
    diff.marker$cluster <- factor(diff.marker$cluster, levels = cluster.order)
  }
  
  # 获取背景颜色数据
  back.data <- map_df(unique(diff.marker$cluster), function(x) {
    tmp <- diff.marker %>%
      filter(cluster == x)
    
    new.tmp <- data.frame(
      cluster = x,
      min = min(tmp$avg_log2FC) - 0.2,
      max = max(tmp$avg_log2FC) + 0.2
    )
    return(new.tmp)
  })
  
  # 获取top genes
  top.marker.tmp <- diff.marker %>%
    group_by(cluster)
  
  top.marker.max <- top.marker.tmp %>%
    slice_max(n = topGeneN, order_by = avg_log2FC)
  
  top.marker.min <- top.marker.tmp %>%
    slice_min(n = topGeneN, order_by = avg_log2FC)
  
  top.marker <- rbind(top.marker.max, top.marker.min)
  
  # 绘图
  p <- ggplot(diff.marker, aes(x = cluster, y = avg_log2FC)) +
    geom_jitter(aes(color = type), size = 0.75) +
    scale_color_manual(values = c("sigDown" = aesCol[1], "sigUp" = aesCol[2])) +
    geom_tile(data = back.data, aes(x = cluster, y = 0, fill = cluster),
              color = "black", height = log2FC.cutoff * 2, alpha = 0.3, show.legend = FALSE) +
    scale_fill_manual(values = tile.col) +
    ggrepel::geom_text_repel(
      data = top.marker,
      aes(x = cluster, y = avg_log2FC, label = gene),
      max.overlaps = 50, ...
    ) +
    theme_classic(base_size = base_size) +
    theme(
      panel.grid = element_blank(),
      legend.position = legend.position,
      legend.title = element_blank(),
      legend.background = element_blank()
    ) +
    xlab("Clusters") + ylab("Average log2FoldChange") +
    guides(color = guide_legend(override.aes = list(size = 5)))
  
  return(p)
}

# 使用调整后的函数绘制火山图
combined_volcano_plot <- jjVolcano_combined()
print(combined_volcano_plot)

