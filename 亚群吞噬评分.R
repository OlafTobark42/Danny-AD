# 加载必要的库
library(Seurat)
library(GSEABase)
library(ggplot2)
# 安装 msigdbr 包（如果尚未安装）
if (!requireNamespace("msigdbr", quietly = TRUE)) {
  install.packages("msigdbr")
}

# 安装 msigdbr 包（如果尚未安装）
if (!requireNamespace("msigdbr", quietly = TRUE)) {
  install.packages("msigdbr")
}

# 加载 msigdbr 包
library(msigdbr)
library(dplyr)

# 从 msigdbr 中获取小鼠的吞噬作用相关基因集
phagocytosis_genes <- msigdbr(species = "Mus musculus", category = "C5", subcategory = "GO:BP") %>%
  filter(grepl("phagocytosis", gs_name, ignore.case = TRUE)) %>%
  pull(gene_symbol)

# 确认获取的基因集
print(phagocytosis_genes)

# 2. 使用 AddModuleScore 计算评分
seu <- AddModuleScore(seu, features = list(phagocytosis_genes), name = "PhagocytosisScore")

# 3. 绘制 VlnPlot
VlnPlot(seu, features = "PhagocytosisScore1", group.by = "group", split.by = "RNA_snn_res.0.2") +
  theme_minimal() +
  labs(title = "Phagocytosis Score by Cluster", x = "Cluster", y = "Phagocytosis Score")


# 多配色方案 -------------------------------------------------------------------
# 加载必要的库
library(Seurat)
library(ggplot2)
library(ggsci)
library(RColorBrewer)

# 创建保存结果的文件夹
output_dir <- "visual/吞噬评分"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 设置需要使用的配色方案，指定来源和对应的方案
color_schemes <- list(
  ggsci = c("npg", "aaas", "lancet", "jama", "nejm", "ucscgb", "d3", "tron", "futurama", "rickandmorty"),
  RColorBrewer = c("Set1", "Set2", "Set3", "Paired", "Dark2", "Accent", "Spectral", "RdYlBu", "RdYlGn", "BrBG")
)

# 评分的列名
score_column <- "PhagocytosisScore1"

# 循环生成和保存图表
for (source in names(color_schemes)) {
  for (scheme in color_schemes[[source]]) {
    
    # 绘制小提琴图
    plot <- VlnPlot(seu, features = score_column, pt.size = 0, group.by = "RNA_snn_res.0.2") +
      stat_summary(fun = median, geom = "crossbar", width = 0.5, color = "black", linewidth = 0.3) +
      theme_minimal() +
      labs(title = paste("Phagocytosis Score -", scheme), x = "Cluster", y = "Phagocytosis Score") +
      theme(legend.position = "none")+
      theme_minimal(base_line_size = 0)
    
    # 动态应用不同的配色
    if (source == "ggsci") {
      plot <- plot + do.call(paste0("scale_fill_", scheme), list())
    } else if (source == "RColorBrewer") {
      plot <- plot + scale_fill_brewer(palette = scheme)
    }
    
    # 保存为 PDF 文件
    ggsave(filename = file.path(output_dir, paste0("PhagocytosisScore_", scheme, ".pdf")), 
           plot = plot, width = 8, height = 6)
  }
}

print(paste("All plots saved to:", output_dir))

