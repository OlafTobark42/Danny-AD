library(ComplexHeatmap)
library(circlize)
library(viridis)

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

# 展平marker基因并保持顺序
all_marker_genes <- unlist(marker_genes)

# 保留在数据中存在的基因
all_marker_genes <- all_marker_genes[all_marker_genes %in% rownames(microglia_cells)]

# 提取对应亚群标签
gene_subgroup <- rep(names(marker_genes), sapply(marker_genes, length))
gene_subgroup <- gene_subgroup[all_marker_genes %in% rownames(microglia_cells)]

# 提取表达数据
expr_matrix <- GetAssayData(microglia_cells, slot = "data")[all_marker_genes, ]

# 计算每个基因在各cluster中的平均表达值
avg_exp <- AverageExpression(microglia_cells, features = all_marker_genes, group.by = "seurat_clusters")
avg_exp_mat <- as.matrix(avg_exp$RNA)

# 对表达矩阵进行 log1p 变换，以减少极端高值的影响
avg_exp_mat <- log1p(avg_exp_mat)

# 确保 avg_exp_mat 的行顺序与 selected_genes 一致
avg_exp_mat <- avg_exp_mat[as.character(all_marker_genes), ]

# 对每个基因（行）进行 Z-score 标准化
avg_exp_mat <- t(scale(t(avg_exp_mat)))

# 制作右侧注释
row_anno <- rowAnnotation(
  Subtype = gene_subgroup,
  col = list(Subtype = structure(viridis(length(unique(gene_subgroup))),
                                 names = unique(gene_subgroup)))
)



# 更新 gene_subgroup，确保顺序和长度与 all_marker_genes 一致
gene_subgroup <- rep(names(marker_genes), sapply(marker_genes, length))
gene_subgroup <- gene_subgroup[unlist(marker_genes) %in% all_marker_genes]
# gene_subgroup <- factor(gene_subgroup, levels = c("HM", "IFAM","IRM","ARM","PAM","MYEM","DAM"))


# 转换 gene_subgroup 为因子，按字母顺序排序
gene_subgroup <- factor(gene_subgroup, levels = sort(unique(gene_subgroup)))

# 确保 gene_subgroup 和 avg_exp_mat 的行数一致
if (length(gene_subgroup) != nrow(avg_exp_mat)) {
  stop("Error: Length of gene_subgroup does not match rows in avg_exp_mat.")
}

# 绘制热图
Heatmap(
  avg_exp_mat[all_marker_genes, ],  # 只取匹配的基因
  name = "Log-scaled\nZ-score",
  col = viridis(100),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_split = gene_subgroup,  # 按亚群分块
  row_title = NULL,
  right_annotation = rowAnnotation(
    Subtype = gene_subgroup,
    col = list(Subtype = structure(
      viridis(length(unique(gene_subgroup))),
      names = unique(gene_subgroup)
    ))
  ),
  column_names_side = "bottom",
  row_names_side = "right",
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 10),
  heatmap_legend_param = list(
    title = "Log-scaled\nZ-score",
    at = c(-2, 0, 2)
  ),
  cell_fun = function(j, i, x, y, w, h, col) {
    grid.rect(x, y, unit(1, "mm"), unit(1, "mm"), gp = gpar(fill = col, col = NA))
  }
)
