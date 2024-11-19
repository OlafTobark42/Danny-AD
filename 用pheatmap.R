library(pheatmap)
library(viridis)
library(Seurat)

# 定义Marker基因
marker_genes <- list(
  HM = c("Cx3cr1", "Tmem119", "P2ry12", "Csf1r", "Zfhx3", "maf", "Cacnb2"),
  IFAM = c("Ccl4", "Ccl3", "Il1b", "Zfp36", "Lgals3", "Cst7", "Atf3", "Jun", "Tnf", "Il1a"),
  IRM = c("Ifitm3", "Rtp4", "Oasl2", "Stat1", "Stat4", "Ifit2", "Irf7", "Oas1", "Plscr1", "Ifit3", "Ifi204"),
  ARM = c("H2-Aa", "H2-Ab1", "H2-Eb1", "Cd74", "Ciita", "H2-D1", "H2-K1", "Nlrc5"),
  PAM = c("Clec7a", "Spp1", "Gpnmb", "Igf1", "Igtax", "Fabp5", "Lgals3"),
  MYEM = c("Lpl", "Cst7", "Igf1", "Spp1", "Fabp5", "Itgax"),
  DAM = c("Spp1", "Gpnmb", "Igf1", "Clec7a", "Lpl", "Lgals3", "Itgax", "Apoe", "Tyrobp", "Cd9", "Cst7", "Cd63", "Fabp5")
)

# 展平 marker_genes 并保持顺序
all_marker_genes <- unlist(marker_genes)

# 提取表达矩阵并筛选 marker_genes
all_marker_genes <- all_marker_genes[all_marker_genes %in% rownames(microglia_cells)]
expr_matrix <- GetAssayData(microglia_cells, slot = "data")[all_marker_genes, ]

# 计算每个基因在各 cluster 的平均表达值
avg_exp <- AverageExpression(microglia_cells, features = all_marker_genes, group.by = "seurat_clusters")
avg_exp_mat <- as.matrix(avg_exp$RNA)

# 对表达矩阵进行 log1p 变换
avg_exp_mat <- log1p(avg_exp_mat)

# 确保 avg_exp_mat 的行顺序与 selected_genes 一致
avg_exp_mat <- avg_exp_mat[as.character(all_marker_genes), ]

# 对每个基因进行 Z-score 标准化
avg_exp_mat <- t(scale(t(avg_exp_mat)))



# 对 avg_exp_mat 的行名进行唯一化
rownames(avg_exp_mat) <- make.unique(rownames(avg_exp_mat))

# 重新生成 gene_subgroup，确保匹配行名
gene_subgroup <- rep(names(marker_genes), sapply(marker_genes, length))
gene_subgroup <- gene_subgroup[unlist(marker_genes) %in% rownames(avg_exp_mat)]

# 创建亚群注释表
annotation_row <- data.frame(Subtype = factor(gene_subgroup, levels = unique(gene_subgroup)))
rownames(annotation_row) <- rownames(avg_exp_mat)  # 此时行名已唯一化
annotation_row$Subtype
# 创建颜色映射
annotation_colors <- list(
  Subtype = c(
    HM = "#E41A1C",     # 红色
    IFAM = "#377EB8",   # 蓝色
    IRM = "#4DAF4A",    # 绿色
    ARM = "#984EA3",    # 紫色
    PAM = "#FF7F00",    # 橙色
    MYEM = "#FFFF33",   # 黄色
    DAM = "#A65628"     # 棕色
  )
)

# 检查结果
print(annotation_colors)


# 绘制热图
pheatmap(
  avg_exp_mat,
  color = viridis(100),
  cluster_rows = FALSE,  # 不对基因聚类
  cluster_cols = FALSE,  # 不对列聚类
  annotation_row = annotation_row,  # 添加亚群注释
  annotation_colors = annotation_colors,
  show_rownames = TRUE,
  show_colnames = TRUE,
  cellwidth = 6,  # 方块宽度
  cellheight = 6,  # 方块高度
  fontsize_row = 6,  # 基因名字体大小
  fontsize_col = 6, # Cluster名字体大小
  main = "Expression of Selected Genes in Microglia Subtypes"
)
